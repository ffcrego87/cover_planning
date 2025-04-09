function [MSTMap, Gr, T, G] = generateMSTMap(P, L, D, ani, d, Lsize, Hsize, nP, npol, MagFactor, debug, MSTHndl)
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
%   debug     - Flag booleana para mostrar figura
%   MSTHndl   - Handle da figura de debug para o mapa MST
%
% Saídas:
%   MSTMap    - Imagem RGB com a árvore geradora desenhada
%   Gr        - Grafo completo com todos os nós válidos
%   T         - Árvore geradora mínima resultante (grafo)
%   G         - Coordenadas dos nós válidos (após remoção nos obstáculos)

map='hex';
x = (d/2):(2*d):Hsize; nx = length(x);
h = 2*d/3;
y = (d/2):(d/sqrt(3)):Lsize; ny = length(y);
[Y,X] = meshgrid(x,y);% A revoir (les notations, x vs y)
Y(1:2:end,:)=Y(1:2:end,:)+d;
% Efface les noeuds dans les polygones
Z = true(ny,nx);
for i=1:npol
    Pi1 = P(L==i,:);
    Pi = [Pi1;Pi1(1,:)];
    for j=1:nx
        for k=1:ny
            if inpolygon(X(k,j),Y(k,j),Pi(:,1),Pi(:,2))
                Z(k,j) = false;
            end
        end
    end
end
x = reshape(X,nx*ny,1);
y = reshape(Y,nx*ny,1);
z = reshape(Z,nx*ny,1);
G = [x(z),y(z)];
% Calcule la matrice d'ajacence
nn = length(G);
Adj = zeros(nn,nn);
%rng('shuffle');
penalty = 1+ani/100;
for j=1:nn
    if debug
        plot(G(j,1),G(j,2),'or')
    end
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
            Adj(j,k) = Lsize+Hsize+sqrt((G(j,1)-G(k,1))^2 + penalty*(G(j,2)-G(k,2))^2)+(1-ani)*randn(1,1)/100;
        end
        Adj(k,j) = Adj(j,k);
    end
end
Gr = graph(Adj);
% Arbre couvrant minimal
T = minspantree(Gr);
%
if debug
    for j=1:height(T.Edges)
        plot([G(T.Edges.EndNodes(j,1),1),G(T.Edges.EndNodes(j,2),1)],[G(T.Edges.EndNodes(j,1),2),G(T.Edges.EndNodes(j,2),2)],'k-')
    end
end
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
%         % Il faut faire quelque chiose mais "l'astuce du siècle" szemble
%         % résoudre le problème, au moins dans certains cas.
%         % Sinon, j'avais imaginé de joindre les composantes connexes...
%         [bins,binsizes] = conncomp(T);
%         newConnection = [0,0,Inf];
%         for j=1:height(Gr.Edges)
%             if bins(Gr.Edges.EndNodes(j,1)) ~= bins(Gr.Edges.EndNodes(j,2)) % not in the same connected component
%                 w=sqrt((G(j,1)-G(k,1))^2 + 1.01*(G(j,2)-G(k,2))^2);
%                 if w<newConnection(3)
%                     newConnection = [Gr.Edges.EndNodes(j,:),Gr.Edges.Weight(j)];
%                 end
%             end
%         end
%         T = addedge(T, newConnection(1), newConnection(2), newConnection(3));
%         fprintf('add : %d - %d , W=%5.2f... ',newConnection(1), newConnection(2), newConnection(3))
    else
        break
    end
end
%%
%    for j=1:height(T.Edges)
%        plot([G(T.Edges.EndNodes(j,1),1),G(T.Edges.EndNodes(j,2),1)],[G(T.Edges.EndNodes(j,1),2),G(T.Edges.EndNodes(j,2),2)],'k-')
%    end
treePath = zeros(height(T.Edges),4);
for j=1:height(T.Edges)
    treePath(j,:) = [G(T.Edges.EndNodes(j,1),1), G(T.Edges.EndNodes(j,1),2), G(T.Edges.EndNodes(j,2),1), G(T.Edges.EndNodes(j,2),2)];
end
MSTMap = 255*ones(Hsize*MagFactor,Lsize*MagFactor,3,'uint8');
MSTMap = insertShape(MSTMap,'Line',treePath*MagFactor,'Color','black','opacity',1);  
if debug
    figure(MSTHndl);
    imagesc(MSTMap)
    axis xy, axis off, axis equal
end
end