function [pgons, npol] = recoverPolygonsFromPL(P, L)
    % npol é o maior rótulo em L
    npol = max(L);
    
    % Prealoca a cell array para armazenar cada polígono
    pgons = cell(npol,1);
    
    % Loop pelos polígonos
    for j = 1:npol
        % Localiza os pontos cujo rótulo seja j
        idx = find(L == j);
        
        % Extrai as coordenadas (x, y) desses pontos
        polyPoints = P(idx, :);
        
        % Se quiser garantir que o primeiro ponto se repita no final
        % (fechando o polígono), faça:
        % polyPoints = [polyPoints; polyPoints(1,:)];
        
        % Armazena no cell array; você pode armazenar direto como Nx2,
        % ou usar 'polyshape' do MATLAB:
        % pgons{j} = polyshape(polyPoints(:,1), polyPoints(:,2));
        pgons{j} = polyPoints;
    end
end
