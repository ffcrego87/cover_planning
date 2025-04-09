function [ALOP, UAR, Surface] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, radius, path, pgons, npol, t, pdf_filename)
    % Inicializar imagem branca
    UARmap = 255 * ones(Hsize * MagFactor, Lsize * MagFactor, 3, 'uint8');
    
    % Desenhar caminho principal
    UARmap = insertShape(UARmap, 'Polygon', path(:)' * MagFactor, 'Color', 'black', 'Opacity', 1);

    % Converter para binário
    BW4 = 255 * uint8(rgb2gray(UARmap) < 128);

    % Dilatação com máscara circular manual
    scale = MagFactor;
    dil_radius = round(scale * radius / 4);
    [x, y] = meshgrid(-dil_radius:dil_radius, -dil_radius:dil_radius);
    circ_mask = (x.^2 + y.^2) <= dil_radius^2;
    SE1 = strel(circ_mask);  % elemento estruturante circular verdadeiro

    BW4 = imdilate(BW4, SE1) / 2;

    % Preencher polígonos
    for j = 1:npol
        polyj = pgons{j}';
        BW4 = insertShape(BW4, 'FilledPolygon', MagFactor * polyj(:)', 'Color', 'white', 'Opacity', 1);
    end

    % Converter novamente para binário
    BW4 = rgb2gray(BW4);

    % Cálculo de métricas
    Uncov  = sum(BW4(:) == 0);
    Island = sum(BW4(:) == 255);
    Surface = length(BW4(:)) - Island;
    UAR = 100 * (Uncov / Surface);
    Surface = Surface / (MagFactor^2);
    ALOP = t(end) * radius / Surface;

    fprintf('Uncovered area ratio (UAR)  = %5.2f%%\n', UAR);
    fprintf('Area length over surface (ALOP) = %.4f\n', ALOP);

    % Plot
    figure('Color','w')
    imagesc(BW4), axis xy off equal
    colormap(gray)
    hold on
    plot(path(1,:) * MagFactor, path(2,:) * MagFactor, 'r-', 'LineWidth', 1.2)

    % Guardar como PDF
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'Auto', ...
             'PaperUnits', 'Inches', ...
             'PaperSize', [pos(3), pos(4)]);
    print(gcf, pdf_filename, '-dpdf', '-r0');
    close(gcf)
end
