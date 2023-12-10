clc
clear all
N=492;debug=true;
% N = 491;debug=true; % spanning tree fail ?
rng(N)
%% Generate polygons write a func   
figure(1)
set(1,'name','main figure : map and path')
clf
title(sprintf('Map number %d',N))
axis off
hold on
% Size of the map
ct=sqrt(3)/2;
d = 5;
Hsize=9*ct*d;
Lsize=9*d;
plot([0,Lsize,Lsize,0,0],...
     [0,0,Hsize,Hsize,0],'k')
sigma = 10;% Rayon moyen des obstacles
nsom = 12;% Nombre maximum de sommets
P = zeros(1000,2);% Nodes of polygons
%P(1:4,:) = [0,0;Hsize,0;Hsize,Lsize;0,Lsize];
nP = 0;% Nombre de points dans P
D = zeros(1000,2);% Edges
%D(1:4,:) = [1:nP;2:nP,1]';
L = zeros(1000,1);% Liste des labels des polygones
nk = 0;
npol = 3;
for j=1:npol
    q = 0;
    nQ = 1;
    while q<=nQ
        while q == 0
            nk = nk+1;
            x=Lsize*rand();
            y=Hsize*rand();
            rho = sigma*rand(nsom,1);
            theta = 2*pi*rand(nsom,1);
            Q = [x+rho.*cos(theta),y+rho.*sin(theta)];
            Q0 = min(Q);
            Q1 = max(Q);
            if min(Q0)>0 && Q1(1)<Lsize && Q1(2)<Hsize
                [k,av] = convhull(Q);
                Q = Q(k,:);
                nQ = length(Q)-1;
                q = 1;
                p = 1;
                % Test if two polygons are included
                i=1;
                included = false;
                while i<j && not(included)
                    Pi1 = P(L==i,:);
                    Pi = [Pi1;Pi1(1,:)];
                    included = inpolygon(Q(1,1),Q(1,2),Pi(:,1),Pi(:,2));
                    included = included || inpolygon(Pi(1,1),Pi(1,2),Q(:,1),Q(:,2));
                    i = i+1;
                end
                if included
                    q = 0;
                end
            else
                q=0;
            end
        end
        if j>1
            q1 = mod(q,nQ)+1;
            R = GetIntersection(P,D(p,1),D(p,2),Q,q,q1);
            if isempty(R)
                p = p+1;
                if p>nP
                    p = 1;
                    q = q+1;
                end
            else
                q = 0;
            end
        else
            p = 1;
            q = q+1;
        end
    end
    P(nP+(1:nQ),:) = Q(1:nQ,:);% Only nodes
    D(nP+(1:nQ),:) = [nP+(1:nQ) ; nP+(2:nQ), nP+1 ]';
    L(nP+(1:nQ)) = j;
    nP = nP+nQ;
    pgon = plot(polyshape(Q(:,1),Q(:,2)));
end
axis([0,Lsize,0,Hsize])
P=P(1:nP,:);
D=D(1:nP,:);
L=L(1:nP);
if debug
    fprintf('Nombre de polygones générés %d\n',nk)
end
%%
x = d:2*d:(Lsize-2*d); nx = length(x);
x(1:2:end)=x(1:2:end);
ct=sqrt(3)/2;%try withiut .5
y = ct*d/2:ct*d:Hsize-ct*d/2; ny = length(y);
[X,Y] = meshgrid(x,y);
X(1:2:end,:)=X(1:2:end,:)+d;
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
            Adj(j,k) = sqrt((G(j,1)-G(k,1))^2 + (G(j,2)-G(k,2))^2);
        else
            Adj(j,k) = Inf;
        end
        Adj(k,j) = Adj(j,k);
    end
end
Gr = graph(Adj);
% Arbre couvrant minimal
T = minspantree(Gr);
if debug
    for j=1:height(T.Edges)
        plot([G(T.Edges.EndNodes(j,1),1),G(T.Edges.EndNodes(j,2),1)],[G(T.Edges.EndNodes(j,1),2),G(T.Edges.EndNodes(j,2),2)],'k-')
    end
end
%%
% hold on
% for j=1:height(T.Edges)
%     plot([G(T.Edges.EndNodes(j,1),1),G(T.Edges.EndNodes(j,2),1)],[G(T.Edges.EndNodes(j,1),2),G(T.Edges.EndNodes(j,2),2)],'k-','LineWidth',2)
% end

fhndl = figure(2);
if debug
    set(2,'visible','on');
    set(2,'name','Minimal covering tree')
else
    set(2,'visible','off');
end
axis off
hold on
for j=1:height(T.Edges)
    plot([G(T.Edges.EndNodes(j,1),1),G(T.Edges.EndNodes(j,2),1)],[G(T.Edges.EndNodes(j,1),2),G(T.Edges.EndNodes(j,2),2)],'k-','LineWidth',2)
end
xlim([0,Lsize])
ylim([0,Hsize])
%%
axhndl = get(fhndl,'Children');
F = getframe(axhndl);
I = frame2im(F);
if not(debug)
    close(fhndl)
end
I = rgb2gray(I);
BW = I<128;
%%
scale = 10;
BW1 = imresize(BW,scale*[Hsize Lsize]);
SE1 = strel("disk",round(scale*d/4));
BW2 = imdilate(BW1,SE1);
if debug
    figure(3)
    set(3,'name','Morphological dilatation of Figure 2')
    imshow(BW2)
end
%% Add obstacles
Img = 255*uint8(BW2);
%Obs = 0*uint8(BW2);
for i=1:npol
    Pi = scale*P(L==i,:);
    Pi(:,2) = scale*Hsize-Pi(:,2);
    Pi=Pi';
    Img = insertShape(Img,'FilledPolygon',{Pi(:)'},...
    'Color', {'black'},'Opacity',1,'LineWidth',1);
%    Obs = insertShape(Obs,'FilledPolygon',{Pi(:)'},...
%    'Color', {'white'},'Opacity',1,'LineWidth',1);
end
Img = rgb2gray(Img);% because of insertShape
%Obs = im2gray(Obs)>128;% because of insertShape
%SE2 = strel("disk",round(1));
%Obs = imdilate(Obs,SE1);
%Img(Obs) = 0;
%%
BW3 = bwmorph(Img,'remove');
BW3 = bwskel(BW3,'MinBranchLength',10);% Almost unnecessary except for some few cases, such as Map N°252 or 221 (requiring MinBranchLength>=8)
if debug
    figure(4)
    set(4,'name','Border of Figure 3')
    imshow(BW3)
end
%% Find parametrized path
Mpath = BW3;
ymax = scale*Hsize;
npath = sum(sum(Mpath));
path = zeros(npath,2);
[yp,xp] = find(Mpath,1);
path(1,:) = [xp,ymax-yp]/scale;
Mpath(yp,xp) = false;
for j=2:npath
    [dy,dx] = find(Mpath(yp-1:yp+1,xp-1:xp+1));
    if isempty(dx)
        if npath-j>0 % strange
             warning('Forgot %d points, probably a matter of connexity of overlapping\nTry increase MinBranchLength in bwskel\n',npath-j)
        end
        break
    end
    dy = dy-2;
    dx = dx-2;
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
%% Plot path on original figure
%use p,D,l and size, d to get that path
%input attirubutes an image out this path down 
figure(1)
plot(path(:,1),path(:,2),'k.')
hold off
fprintf('Length if the path %d\n',2*(height(Gr.Nodes)+1)*d)
%% Explore the tree
% deb =  T.Edges.EndNodes(1,1);
% DFT = dfsearch(T, deb);
% e = DFT(1);
% for j=2:length(DFT)
%     nexte = DFT(j);
%     plot([G(e,1),G(nexte,1)],[G(e,2),G(nexte,2)],'b-','LineWidth',2)
%     e = nexte;
% end
% plot([G(e,1),G(deb,1)],[G(e,2),G(deb,2)],'b-','LineWidth',2)
%%
% t = 0;
% e = T.Edges.EndNodes(1,1);
% for j=1:15
%     text(G(e,1),G(e,2),num2str(e))
%     plot(G(e,1)+d/3*cos(t-pi/2),G(e,2)+d/3*sin(t-pi/2),'bx')
%     nexta = Inf;
%     nb = neighbors(T,e);
%     for k=1:length(nb)
%         f = nb(k);
%         a = atan2(G(f,2)-G(e,2),G(f,1)-G(e,1))-t;
%         if a<0
%             a = a+2*pi;
%         end
%         if a<nexta
%             nexta = a;
%             nexte = f;
%         end
%     end
%     plot([G(e,1),G(nexte,1)]+d/3*cos(t-pi/2),[G(e,2),G(nexte,2)]+d/3*sin(t-pi/2),'b-','LineWidth',2)
%     fprintf('%3d/ %3d ---> %3d (%f ---> %f)\n',j,e,nexte,t,t+a)
%     e = nexte;
%     t = a;
% end
%%
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