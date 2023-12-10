
% Generate random points inside a bounding box
[X,Y]=meshgrid(0:1:10,0:1:10);
x=reshape(X,[],1);
y=reshape(Y,[],1);
% Define polygonal obstacles
obstacles = {[1 1; 1 4; 4 4; 4 1], [6 6; 6 7; 7 7; 7 6]};

% Plot points and obstacles
figure;
hold on;
plot(x,y,'.','MarkerSize',20);
for i = 1:length(obstacles)
    obstacle = obstacles{i};
    patch(obstacle(:,1), obstacle(:,2), 'k');
end
xlim([0 10]);
ylim([0 10]);
axis equal;

% Remove Voronoi pts inside obstacles
for i = 1:length(obstacles)
    obstacle = obstacles{i};
    [in,on] = inpolygon(x,y,obstacle(:,1),obstacle(:,2));
    is_inside = in | on;
    x=x(~is_inside);
    y=y(~is_inside);
    x=[x ; obstacle(:,1)];
    y=[y ; obstacle(:,2)];
end
p=[x y];

% Generate Voronoi diagram
[v,c] = voronoin(p);

% Remove Voronoi vertices outside bounding box
is_outside = v(:,1)<0 | v(:,1)>10 | v(:,2)<0 | v(:,2)>10;
% v(:,1)(is_outside) = NaN;
% v(:,2)(is_outside) = NaN;
index=find(is_outside);
%while(~isempty(index))
for j=1:length(index)
    [v,c] = voronoielimnation(index(j),v,c);
    %index(1)=[];
    index=index-1;
end

% Remove Voronoi edges inside obstacles
for i = 1:length(obstacles)
    obstacle = obstacles{i};
    [in,on] = inpolygon(v(:,1),v(:,2),obstacle(:,1),obstacle(:,2));
    is_inside = in | on;
    index=find(is_inside);
    for j=1:length(index)
        [v,c] = voronoielimnation(index(j),v,c);
        index=index-1;
    end
end

% Plot Voronoi diagram
for j=1:length(c)
    plot([v(c{j},1)' v(c{j}(1),1)],[v(c{j},2)' v(c{j}(1),2)],'r');
end

edges=[];
for z=1:length(v)
    neighbor=[];
    for j=1:length(c)
        if ismember(z,c{j})
%             idx_neighbor=mod(find(z==c{j})+[-2 0],length(c{j}))+1;
%             neighborsj=c{j}(idx_neighbor);
            neighborsj=c{j}(c{j}>z);
            neighbor=[neighbor setdiff(neighborsj,neighbor)];
        end
    end
    edges=[edges;[z*ones(size(neighbor));neighbor]'];
end

idxs = edges;


%%

% Remove non-existent 

num_v = size(v,1);
idxs = nchoosek(1:num_v,2);

dist = hypot(v(idxs(:,1),1) - v(idxs(:,2),1), ...
             v(idxs(:,1),2) - v(idxs(:,2),2));
lendist = length(dist);

for i = 1:length(obstacles)
    obstacle = obstacles{i};
    for j=1:length(idxs)
        [xi,yi] = polyxpoly(v(idxs(j,:),1),v(idxs(j,:),2),obstacle(:,1),obstacle(:,2));
        if ~isempty(xi)
            dist(j)=10^6;
        end
    end
end

G = graph(idxs(:,1),idxs(:,2));

figure
hGraph = plot(G,'XData',v(:,1),'YData',v(:,2),'LineStyle','none','NodeLabel',{});
hold on
% Draw the outside border
%plot(x,y,'r-')
hold off

%%

Aeq = spalloc(num_v,length(idxs),num_v*(num_v-1)); % Allocate a sparse matrix
for ii = 1:num_v
    whichIdxs = (idxs == ii); % Find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); % Include trips where ii is at either end
    Aeq(ii,:) = whichIdxs'; % Include in the constraint matrix
end
beq = 2*ones(num_v,1);

intcon = 1:lendist;
lb = zeros(lendist,1);
ub = inf(lendist,1);

%opts = optimoptions('intlinprog','Display','off');
[x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub);

x_tsp = logical(round(x_tsp));
Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2));

hold on
highlight(hGraph,Gsol,'LineStyle','-')
title('Solution with Subtours')

%%

tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % number of subtours
fprintf('# of subtours: %d\n',numtours);

A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
b = [];
while numtours > 1 % Repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,lendist,num_v)]; % A guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1) + 1; % Counter for indexing
        subTourIdx = find(tourIdxs == ii); % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:size(variations,1)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx) - 1; % One less trip than subtour stops
    end

    % Try to optimize again
    [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,A,b,Aeq,beq,lb,ub);
    x_tsp = logical(round(x_tsp));
    Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2));
    
    % Visualize result
    hGraph.LineStyle = 'none'; % Remove the previous highlighted path
    highlight(hGraph,Gsol,'LineStyle','-')
    drawnow
    
    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % number of subtours
    fprintf('# of subtours: %d\n',numtours)
end

for i = 1:length(obstacles)
    obstacle = obstacles{i};
    patch(obstacle(:,1), obstacle(:,2), 'k');
end

title('Solution with Subtours Eliminated');
hold off

