function [pos1sol,pos2sol,t] = TSP_plan2(Lsize,Hsize,d,P,~,L,speed)

% Generate random points inside a bounding box
% Define polygonal obstacles
obstacles=cell(max(L),1);
for i=1:max(L)
    obstacles{i} = P(L==i,:);
end

%%
[X,Y]=meshgrid(0:(d/2):Lsize,0:(d/2):Hsize);
v=[reshape(X,[],1) reshape(Y,[],1)];

for i = 1:length(obstacles)
    obstacle = obstacles{i};
    [in,on] = inpolygon(v(:,1),v(:,2),obstacle(:,1),obstacle(:,2));
    is_inside = in | on;
    v=v(~is_inside,:);
end

num_v = size(v,1);
idxs = nchoosek(1:num_v,2);

dist = hypot(v(idxs(:,1),1) - v(idxs(:,2),1), ...
             v(idxs(:,1),2) - v(idxs(:,2),2));
lendist = length(dist);

for i = 1:length(obstacles)
    obstacle = obstacles{i};
    for j=1:length(idxs)
        [xi,~] = polyxpoly(v(idxs(j,:),1),v(idxs(j,:),2),obstacle(:,1),obstacle(:,2));
        if ~isempty(xi)
            dist(j)=10^6;
        end
    end
end

Dist_MAT = zeros(num_v);

for i=1:size(idxs,1)
    Dist_MAT(idxs(i,1),idxs(i,2))=dist(i);
    Dist_MAT(idxs(i,2),idxs(i,1))=dist(i);
end

[p,~] = tspsearch(Dist_MAT,10);


pos1sol=v(p,1);
pos2sol=v(p,2);
t = zeros(size(pos1sol));
for i=2:length(t)
    t(i)=t(i-1)+norm([pos1sol(i)-pos1sol(i-1);pos2sol(i)-pos2sol(i-1)])/speed;
end
end