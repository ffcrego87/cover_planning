function [pos1sol,pos2sol,t] = MST_plan_out(Lsize,Hsize,d,P,D,L,Pout,speed)

%%
x = 0:d:Lsize; nx = length(x);
y = 0:d:Hsize; ny = length(y);
nP = 0;% Nombre de points dans P
% Efface les noeuds dans les polygones
Z = true(ny,nx);
for i=1:max(L)
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
% Outer polygon
for j=1:nx
    for k=1:ny
        if ~inpolygon(x(j),y(k),Pout(:,1),Pout(:,2))
            Z(k,j) = false;
        end
    end
end

[X,Y] = meshgrid(x,y);
x = reshape(X,nx*ny,1);
y = reshape(Y,nx*ny,1);
z = reshape(Z,nx*ny,1);
G = [x(z),y(z)];
% Calcule la matrice d'ajacence
nn = length(G);
Adj = zeros(nn,nn);
for j=1:nn
    for k=1:j-1
        p = 1;
        R = [];
        while p<=nP && isempty(R)
            R = GetIntersection(G,j,k,P,D(p,1),D(p,2));
            p = p+1;
        end
        [xi,~] = polyxpoly([G(j,1) G(k,1)],[G(j,2) G(k,2)],Pout(:,1),Pout(:,2));
        if isempty(R) && isempty(xi)
            Adj(j,k) = sqrt((G(j,1)-G(k,1))^2 + (G(j,2)-G(k,2))^2);
        else
            Adj(j,k) = Inf;
        end
        Adj(k,j) = Adj(j,k);
    end
end
% Find nodes with only infinite distances
isolatedNodes = all(isinf(Adj+diag(inf*ones(length(Adj),1))), 2);

% Remove isolated nodes from G and Adj
G(isolatedNodes, :) = [];
Adj(isolatedNodes, :) = [];
Adj(:, isolatedNodes) = [];

Gr = graph(Adj);
% Arbre couvrant minimal
T = minspantree(Gr);

%%
fhndl = figure(2);
set(2,'visible','off');
axis off
hold on
for j=1:height(T.Edges)
    plot([G(T.Edges.EndNodes(j,1),1),G(T.Edges.EndNodes(j,2),1)],[G(T.Edges.EndNodes(j,1),2),G(T.Edges.EndNodes(j,2),2)],'k-','LineWidth',2)
end
xlim([0 Lsize]);
ylim([0 Hsize]);
axhndl = get(fhndl,'Children');
F = getframe(axhndl);
I = frame2im(F);
close(fhndl)
I = rgb2gray(I);
BW = I<128;

%%
scale = 10;
BW1 = imresize(BW,scale*[Hsize Lsize]);
SE1 = strel("disk",round(scale*d/8));
BW2 = imdilate(BW1,SE1);

%% Add obstacles
Img = 255*uint8(BW2);

Pout_aux = scale * Pout;
Pout_aux(:, 2) = scale * Hsize - Pout_aux(:, 2);
Pout_aux = Pout_aux';
polygonMask = poly2mask(Pout_aux(1, :), Pout_aux(2, :), size(Img, 1), size(Img, 2));
outsideMask = ~polygonMask;
Img(outsideMask) = 0;

for i=1:max(L)
    Pi = scale*P(L==i,:);
    Pi(:,2) = scale*Hsize-Pi(:,2);
    Pi=Pi';
    Img = insertShape(Img,'FilledPolygon',{Pi(:)'},...
        'Color', {'black'},'Opacity',1,'LineWidth',1);
    Img = rgb2gray(Img);
end

%%
BW3 = bwmorph(Img,'remove');

%% Find parametrized path
Mpath = BW3;
no=10;
ymax = scale*Hsize;
npath = sum(sum(Mpath));
path = zeros(npath,2);
indices = find(Mpath);
randomIndex = indices(randi(length(indices)));
[yp,xp] = ind2sub(size(Mpath), randomIndex);
path(1,:) = [xp,ymax-yp]/scale;
Mpath(yp,xp) = false;
for j=2:npath
    [dy,dx] = find(Mpath(max(yp-no,1):min(yp+no,size(Mpath,1)),max(xp-no,1):min(xp+no,size(Mpath,2))));
    if isempty(dx)
        path = path(1:(j-1),:);
        break
    end
    dy = dy-min(no+1,yp);
    dx = dx-min(no+1,xp);
    if length(dx)>1 % 4-order neighboor first
        [~,k] = sort(dx.^2+dy.^2);
        dx = dx(k(1));
        dy = dy(k(1));
    end
    yp = yp+dy;
    xp = xp+dx;
    path(j,:) = [xp,ymax-yp]/scale;
    Mpath(yp,xp) = false;
end
pos1sol=path(:,1);
pos2sol=path(:,2);
t = zeros(size(pos1sol));
for i=2:length(t)
    t(i)=t(i-1)+norm([pos1sol(i)-pos1sol(i-1);pos2sol(i)-pos2sol(i-1)])/speed;
end
    
end

function R = GetIntersection(P,i1,i2,Q,j1,j2)
b = (P(i1,1)-P(i2,1))*(Q(j2,2)-Q(j1,2))-(P(i1,2)-P(i2,2))*(Q(j2,1)-Q(j1,1));
a = (Q(j1,1)-P(i2,1))*(Q(j2,2)-Q(j1,2))-(Q(j1,2)-P(i2,2))*(Q(j2,1)-Q(j1,1));
t = a / b;
if t>0 && t<1
    R = [t*P(i1,1)+(1-t)*P(i2,1) , t*P(i1,2)+(1-t)*P(i2,2)];
    s = (R(1)-Q(j1,1))*(Q(j2,1)-Q(j1,1))+(R(2)-Q(j1,2))*(Q(j2,2)-Q(j1,2));
    s = s/((Q(j2,1)-Q(j1,1))^2+(Q(j2,2)-Q(j1,2))^2);
    if s>0 && s<1
        return
    end
end
R = [];
end