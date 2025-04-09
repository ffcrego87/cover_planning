function [MSTMap, Gr, T, G] = generateMSTMap(P, L, D, ani, d, Lsize, Hsize, nP, npol, MagFactor)
% generateMSTMap - Gera uma árvore geradora mínima (MST) considerando obstáculos
%
% Entradas:
%   P         - Conjunto de pontos (Nx2) que define obstáculos
%   L         - Etiquetas correspondentes a cada ponto (para agrupamento por polígono)
%   D         - Par ordenado de pontos que define restrições adicionais
%   ani       - Parâmetro de anisotropia (percentagem)
%   d         - Espaçamento da grelha
%   Lsize     - Largura do mapa
%   Hsize     - Altura do mapa
%   nP        - Número de pares em D
%   npol      - Número de polígonos
%   MagFactor - Fator de ampliação da imagem
%
% Saídas:
%   MSTMap    - Imagem RGB com a árvore geradora desenhada
%   Gr        - Grafo completo com todos os nós válidos
%   T         - Árvore geradora mínima resultante (grafo)
%   G         - Coordenadas dos nós válidos (após remoção nos obstáculos)

% Génère un graphe puis un arbre couvrant
x = d/2:d:Lsize-d/2; nx = length(x);
y = d/2:d:Hsize-d/2; ny = length(y);
% Efface les noeuds dans les polygones
Z = true(ny,nx);
for i=1:npol
    Pi1 = P(L==i,:);
    Pi = [Pi1;Pi1(1,:)];
    for j=1:nx
        for k=1:ny
            if inpolygon(x(j),y(k),Pi(:,1),Pi(:,2))
                Z(k,j) = false;
            end
        end
    end
end
[X,Y] = meshgrid(x,y);
x = reshape(X,nx*ny,1);
y = reshape(Y,nx*ny,1);
z = reshape(Z,nx*ny,1);
G = [x(z),y(z)];
% Calcule la matrice d'adjacence
nn = length(G);
Adj = zeros(nn,nn);
penalty = 1+ani/100;
for j=1:nn
    for k=1:j-1
        p = 1;
        R = [];
        while p<=nP && isempty(R)
            R = GetIntersection(G,j,k,P,D(p,1),D(p,2));
            p = p+1;
        end
        if isempty(R)
            Adj(j,k) = sqrt((G(j,1)-G(k,1))^2 + penalty*(G(j,2)-G(k,2))^2)+(1-ani)*randn(1,1)/100;
        else
            Adj(j,k) = Lsize+Hsize+sqrt((G(j,1)-G(k,1))^2 + penalty*(G(j,2)-G(k,2))^2)+(1-ani)*randn(1,1)/100;% L'astuce du siècle !!!
        end
        Adj(k,j) = Adj(j,k);
    end
end
Gr = graph(Adj);
% Arbre couvrant minimal
T = minspantree(Gr);
% Recherche de non connexité
while 1
    j = 0;
    while j<height(T.Edges)
        j = j+1;
        if T.Edges.Weight(j)>2*(Lsize+Hsize)
            fprintf('%d : %d - %d , W=%5.2f... ',j,T.Edges.EndNodes(j,1),T.Edges.EndNodes(j,2),T.Edges.Weight(j))
            break
        end
    end
    if j<height(T.Edges) % Aïe : plusieurs composantes connexes
        fprintf('removed !\n')
        T = rmedge(T,T.Edges.EndNodes(j,1),T.Edges.EndNodes(j,2));
    else
        break
    end
end
%%
treePath = zeros(height(T.Edges),4);
for j=1:height(T.Edges)
    treePath(j,:) = [G(T.Edges.EndNodes(j,1),1), G(T.Edges.EndNodes(j,1),2), G(T.Edges.EndNodes(j,2),1), G(T.Edges.EndNodes(j,2),2)];
end
MSTMap = 255*ones(Hsize*MagFactor,Lsize*MagFactor,3,'uint8');
MSTMap = insertShape(MSTMap,'Line',treePath*MagFactor,'Color','black','opacity',1);
end