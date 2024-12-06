function [pos1sol,pos2sol,t] = TSP_plan2_out(Lsize,Hsize,d,P,~,L,Pout,speed)

% Generate random points inside a bounding box
% Define polygonal obstacles
obstacles=cell(max(L),1);
for i=1:max(L)
    obstacles{i} = P(L==i,:);
end

%%
[X,Y]=meshgrid(0:(d/2):Lsize,0:(d/2):Hsize);
v=[reshape(X,[],1) reshape(Y,[],1)];

[in,on] = inpolygon(v(:,1),v(:,2),Pout(:,1),Pout(:,2));
is_inside = in | on;
v=v(is_inside,:);

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
for j=1:length(idxs)
    [xi,~] = polyxpoly(v(idxs(j,:),1),v(idxs(j,:),2),Pout(:,1),Pout(:,2));
    if ~isempty(xi)
        path = pathfinder_poly2(v(idxs(j,:),1),v(idxs(j,:),2), Pout);
        diff2 = path(2:end,:)-path(1:(end-1),:);
        dists = sqrt(diff2(:,1).^2+diff2(:,2).^2);
        dist(j)=sum(dists);
    end
end

Dist_MAT = zeros(num_v);

for i=1:size(idxs,1)
    Dist_MAT(idxs(i,1),idxs(i,2))=dist(i);
    Dist_MAT(idxs(i,2),idxs(i,1))=dist(i);
end

[p,custo] = tspsearch(Dist_MAT,10);

disp(custo)

pos1sol=v(p,1);
pos2sol=v(p,2);

% waypoints = [];
% for i=1:(length(pos1sol_prov)-1)
%     path = pathfinder_poly2([pos1sol_prov(i) pos2sol_prov(i)], [pos1sol_prov(i+1) pos2sol_prov(i+1)], Pout);
%     waypoints = [waypoints;path(1:(end-1),:)];
% end
% waypoints = [waypoints;[pos1sol_prov(end) pos2sol_prov(end)]];
% 
% pos1sol = waypoints(:,1)';
% pos2sol = waypoints(:,2)';

t = zeros(size(pos1sol));
for i=2:length(t)
    t(i)=t(i-1)+norm([pos1sol(i)-pos1sol(i-1);pos2sol(i)-pos2sol(i-1)])/speed;
end
end