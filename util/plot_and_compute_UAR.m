function [ALOP, UAR, Surface] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, radius, path, pgons, npol, t, speed, pdf_filename)

    %---------------------------------------
    % 1) Definir parâmetros e criar padding
    %---------------------------------------
    % Fator de escala e raio do disco
    scale = MagFactor;
    dil_radius = round(scale * radius / 4);

    % Quantidade de padding (margem extra)
    padSize = dil_radius + 2;

    % Dimensões da imagem com padding
    paddedH = (Hsize + 2*padSize) * MagFactor;
    paddedL = (Lsize + 2*padSize) * MagFactor;

    % Criar imagem branca (255) com padding
    UARmap_padded = 255 * ones(paddedH, paddedL, 3, 'uint8');

    %---------------------------------------
    % 2) Transladar (shift) o caminho
    %---------------------------------------
    % Ajustar o caminho para dentro da imagem maior
    pathShifted = path + padSize;  % soma padSize a x e y

    % Desenhar caminho na imagem com padding
    UARmap_padded = insertShape(UARmap_padded, ...
                                'Polygon', ...
                                (pathShifted(:)' * MagFactor), ...
                                'Color', 'black', ...
                                'Opacity', 1);

    % Converter para binário
    BW4_padded = 255 * uint8(rgb2gray(UARmap_padded) < 128);

    %---------------------------------------
    % 3) Dilatação com máscara circular real
    %---------------------------------------
    [xx, yy] = meshgrid(-dil_radius:dil_radius, -dil_radius:dil_radius);
    circ_mask = (xx.^2 + yy.^2) <= dil_radius^2;   % máscara circular
    SE1 = strel(circ_mask);                        % strel circular
    BW4_padded = imdilate(BW4_padded, SE1) / 2;    % dilata

    %---------------------------------------
    % 4) Preencher polígonos (também shift)
    %---------------------------------------
    for j = 1:npol
        polyj = pgons{j}';
        polyjShifted = polyj + padSize; % transladar vértices
        BW4_padded = insertShape(BW4_padded, ...
                                 'FilledPolygon', ...
                                 (polyjShifted(:)' * MagFactor), ...
                                 'Color', 'white', ...
                                 'Opacity', 1);
    end

    % Converter para binário
    BW4_padded = rgb2gray(BW4_padded);

    %---------------------------------------
    % 5) "Recortar" (unpad) de volta ao tamanho
    %---------------------------------------
    rowStart = padSize * MagFactor + 1;
    rowEnd   = rowStart + (Hsize * MagFactor) - 1;
    colStart = padSize * MagFactor + 1;
    colEnd   = colStart + (Lsize * MagFactor) - 1;

    BW4 = BW4_padded(rowStart:rowEnd, colStart:colEnd);

    %---------------------------------------
    % 6) Cálculo de métricas
    %---------------------------------------
    Uncov  = sum(BW4(:) == 0);
    Island = sum(BW4(:) == 255);
    Surface = numel(BW4) - Island;
    UAR = 100 * (Uncov / Surface);
    Cov = (Surface-Uncov) / (MagFactor^2);
    ALOP = (t(end)/speed) * radius / Cov;

    fprintf('Uncovered area ratio (UAR)  = %5.2f%%\n', UAR);
    fprintf('Area length over surface (ALOP) = %.4f\n', ALOP);

    %---------------------------------------
    % 7) Plot
    %---------------------------------------
    figure('Color','w')
    imagesc(BW4), axis xy off equal
    colormap(gray)
    hold on

    % Para desenhar o caminho em vermelho na figura final,
    % basta usar o caminho original (sem pad), mas escalado:
    pathPlot = path * MagFactor;
    plot(pathPlot(1,:), pathPlot(2,:), 'r-', 'LineWidth', 1.2)

    %---------------------------------------
    % 8) Guardar como PDF
    %---------------------------------------
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'Auto', ...
             'PaperUnits', 'Inches', ...
             'PaperSize', [pos(3), pos(4)]);
    print(gcf, pdf_filename, '-dpdf', '-r0');
    close(gcf)
end
