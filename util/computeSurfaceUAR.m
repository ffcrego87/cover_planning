function [Surface, UAR] = computeSurfaceUAR(Hsize, Lsize, MagFactor, path, radius, npol, pgons)
    % Cria imagem base (mapa)
    UARmap = 255 * ones(Hsize * MagFactor, Lsize * MagFactor, 3, 'uint8');
    
    % Inserir o polígono inicial
    patht = path';
    UARmap = insertShape(UARmap, 'Polygon', patht(:)' * MagFactor, ...
                         'Color', 'black', 'Opacity', 1);
    
    % Converte para grayscale e binariza
    BW4 = 255 * uint8(rgb2gray(UARmap) < 128);
    scale = MagFactor;
    
    % Dilata a região
    SE1 = strel("disk", round(scale * radius / 4));
    BW4 = imdilate(BW4, SE1) / 2;
    
    % Preenche as regiões dos polígonos adicionais
    for j = 1:npol
        polyj = pgons{j}';
        BW4 = insertShape(BW4, 'FilledPolygon', MagFactor * polyj(:)', ...
                          'Color', 'white', 'Opacity', 1);
    end
    
    % Converte para grayscale (necessário porque insertShape retorna RGB)
    BW4 = rgb2gray(BW4);
    
    % Cálculo de áreas
    Uncov  = sum(BW4(:) == 0);       % Pixels descobertos (pretos)
    Island = sum(BW4(:) == 255);     % Pixels brancos (fora da região)
    
    % Área total (em pixels) menos as áreas fora do polígono
    Surface = numel(BW4) - Island;
    
    % Uncovered Area Ratio
    UAR = 100 * (Uncov / Surface);
    
    % Ajusta a área para compensar o fator de ampliação
    Surface = Surface / (MagFactor^2);
    
    % Exibe o valor no console (opcional)
    fprintf('Uncovered area ratio (UAR)  = %5.2f%%\n', UAR);
end
