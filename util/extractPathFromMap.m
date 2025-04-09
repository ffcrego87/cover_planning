function [path, LOP] = extractPathFromMap(MSTMap, pgons, d, MagFactor, Lsize, Hsize, npol, Gr, N, debug, mainPlot, mainMap, DilatHndl, BorderHndl, mainHndl)
% extractPathFromMap - Extrai um caminho esquelético de uma imagem com obstáculos
%
% Entradas:
%   MSTMap      - Mapa original RGB
%   pgons       - Célula com polígonos de obstáculos (1:npol)
%   d           - Parâmetro de escala base (por ex., largura de aresta)
%   MagFactor   - Fator de ampliação da imagem
%   Lsize, Hsize- Dimensões do mapa original
%   npol        - Número de polígonos (obstáculos)
%   Gr          - Grafo associado à estrutura (para o LOP)
%   N           - Número da instância do mapa (para título)
%   debug       - Flag booleana para ativar figuras de depuração
%   mainPlot    - Flag booleana para gerar imagem final
%   mainMap     - Mapa de fundo onde o caminho será desenhado
%   DilatHndl   - Handle da figura de dilatação
%   BorderHndl  - Handle da figura da fronteira
%   mainHndl    - Handle da figura principal
%
% Saídas:
%   path        - Caminho extraído (Nx2) com coordenadas em escala original
%   LOP         - Length of Path (valor estimado)

I = rgb2gray(MSTMap);
BW = I<128;
BW1 = BW;
scale = MagFactor;
% ALTERNATIVE 2025 : END
SE1 = strel("disk",round(scale*d/4));
BW2 = imdilate(BW1,SE1);
if debug
    figure(DilatHndl);
    imagesc(255-BW2),axis xy, axis off, axis equal
    colormap(gray)
end
% Add obstacles
Img = 255*uint8(BW2);
for j=1:npol
    polyj = pgons{j}';
    Img = insertShape(Img,'FilledPolygon',MagFactor*polyj(:)','Color','black','opacity',1);
end
Img = rgb2gray(Img);% because of insertShape
% Remove interior pixels (keep only border)
BW3 = bwmorph(Img,'remove');
BW3 = bwskel(BW3,'MinBranchLength',20);% Almost unnecessary except for some few cases, such as Map N°252 or 221 (requiring MinBranchLength>=8)
if debug
    figure(BorderHndl);
    imagesc(255-BW3),axis xy, axis off, axis equal
    colormap(gray)
end
% Find parametrized path
Mpath = BW3;
xmax = scale*Lsize;
ymax = scale*Hsize;
npath = sum(sum(Mpath));
path = zeros(npath,2);
[yp,xp] = find(Mpath,1);
path(1,:) = [xp,yp]/scale;
Mpath(yp,xp) = false;
maxsize = 3*MagFactor;
for j=2:npath
    incr = 1;
    [dy,dx] = find(Mpath(max(1,yp-incr):min(ymax,yp+incr),max(1,xp-incr):min(xmax,xp+incr)));    
    while isempty(dx) && incr<maxsize
        incr = incr+1;
        [dy,dx] = find(Mpath(max(1,yp-incr):min(ymax,yp+incr),max(1,xp-incr):min(xmax,xp+incr)));
    end
    if isempty(dx)      % Pas très normal...
        if  npath-j>100  % ... et ça va se payer !
            warning('Forgot %d points during edge detection.\nTry increasing the neighborhood maxsize=%d\n',npath-j,maxsize)
        end
        break
    end
    dy = dy-2;
    dx = dx-2;
    if length(dx)>1 % 4-order neighboor first
        % Not so common, should be avoided, corresponds to branch...
        [~,k] = sort(dx.^2+dy.^2);
        dx = dx(k(1));
        dy = dy(k(1));
    end
    yp = yp+dy;
    xp = xp+dx;
    path(j,:) = [xp,yp]/scale;
    Mpath(yp,xp) = false;
end


LOP = 2*(height(Gr.Nodes)+1)*d;


fprintf('Length of the path (LOP)    = %d\n',LOP)
% Partie graphique
if mainPlot
    figure(mainHndl)
    patht=path';
    mainMap = insertShape(mainMap,'Polygon',patht(:)'*MagFactor,'Color','black','opacity',1);
    imagesc(mainMap)
    title(sprintf('Map number %d',N))
    axis xy, axis off, axis equal
    drawnow
end

end



