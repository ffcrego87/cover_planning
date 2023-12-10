function [pos1sol,pos2sol,t] = TSP_plan(Lsize,Hsize,d,P,D,L,speed)

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

paths = AllPath(adjacency(Gsol), 1, 2);

if length(paths{1})>2
    path = v(paths{1},:);
else
    path = v(paths{2},:);
end

pos1sol=path(:,1);
pos2sol=path(:,2);
t = zeros(size(pos1sol));
for i=2:length(t)
    t(i)=t(i-1)+norm([pos1sol(i)-pos1sol(i-1);pos2sol(i)-pos2sol(i-1)])/speed;
end
end

function p = AllPath(A, s, t)
% Find all paths from node #s to node #t
% INPUTS:
%   A is (n x n) symmetric ajadcent matrix
%   s, t are node number, in (1:n)
% OUTPUT
%   p is M x 1 cell array, each contains array of
%   nodes of the path, (it starts with s ends with t)
%   nodes are visited at most once.
if s == t
    p =  {s};
    return
end
p = {};
As = A(:,s)';
As(s) = 0;
neig = find(As);
if isempty(neig)
    return
end
A(:,s) = [];
A(s,:) = [];
neig = neig-(neig>=s);
t = t-(t>=s); 
for n=neig
    p = [p; AllPath(A,n,t)]; %#ok
end
p = cellfun(@(a) [s, a+(a>=s)], p, 'unif', 0);
end %AllPath
