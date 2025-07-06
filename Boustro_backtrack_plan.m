function [pos1sol,pos2sol,t] = Boustro_backtrack_plan(Lsize,Hsize,d,P,~,L,speed)
% One–pass boustrophedon + nearest-gap back-tracking (straight-line hop)
%
% INPUT  – identical to previous Boustro_plan
% OUTPUT – xPath, yPath = full path;  t = timestamps at constant speed

%-------------------------------- 1. build free grid --------------------------------
step   = d/2;                                   % 50 % lateral overlap
[X,Y]  = meshgrid(0:step:Lsize, 0:step:Hsize);
gridXY = [X(:)  Y(:)];

obs = cell(max(L),1);
for i = 1:max(L),  obs{i} = P(L==i,:);  end

freeMask = true(size(gridXY,1),1);
for i = 1:numel(obs)
    [in,on]   = inpolygon(gridXY(:,1),gridXY(:,2), obs{i}(:,1),obs{i}(:,2));
    freeMask  = freeMask & ~(in|on);
end
freePts  = gridXY(freeMask,:);                  % candidate way-points
visited  = false(size(freePts,1),1);

% Pre-compute row membership
[rowsY,~,rowIdx] = uniquetol(freePts(:,2), step/100,'DataScale',1);

%-------------------------------- core loop -----------------------------------------
dirY    =  1;                                   % 1 = rows ↑,  -1 = rows ↓
path    = zeros(0,2);                           % will be grown
currPos = freePts(1,:);                         % start in first free cell

while any(~visited)
    % -------- sweep pass along rows in current dirY --------
    if dirY>0,  rowOrder = 1:numel(rowsY); else, rowOrder = numel(rowsY):-1:1; end
    for kRow = 1:numel(rowOrder)
        r       = rowOrder(kRow);
        ptsIdx  = find(rowIdx==r & ~visited);   % only still-unvisited pts on row
        if isempty(ptsIdx),  continue,  end
        
        % sort along X and reverse on even kRow to make ←→ zig-zag
        rowPts  = freePts(ptsIdx,:);
        rowPts  = sortrows(rowPts,1,'ascend');
        if mod(kRow,2)==0, rowPts = flipud(rowPts); end
        
        % ---- split by gaps (>1.1*step) to force “turn back” at obstacles
        gapMask     = [true; abs(diff(rowPts(:,1))) > 1.1*step];
        segLabel    = cumsum(gapMask);

        % find closest segment
        [~,ix]   = min( hypot(rowPts(:,1)-currPos(1), rowPts(:,2)-currPos(2)) );
        s   = segLabel(ix);
        
        seg  = rowPts(segLabel==s,:);
        path      = [path ; seg];           %#ok<AGROW>
        visited( ismember(freePts,seg,'rows') ) = true;
        currPos   = seg(end,:);
    end
    if ~any(~visited), break, end               % coverage complete
    
    % -------- straight-line hop to nearest uncovered point --------
    remIdx   = find(~visited);
    [~,ix]   = min( hypot(freePts(remIdx,1)-currPos(1), freePts(remIdx,2)-currPos(2)) );
    tgtIdx   = remIdx(ix);
    tgtPos   = freePts(tgtIdx,:);
    
    % append the hop and mark start point visited (target will be visited in next pass)
    path     = [path ; tgtPos];    %#ok<AGROW>
    visited(tgtIdx) = true;                     % mark as visited immediately
    currPos  = tgtPos;
    
    dirY     = -dirY;                           % flip sweep direction (backwards)
end

%-------------------------------- path ----------------------------------------

pos1sol_prov = path(:,1);
pos2sol_prov = path(:,2);

waypoints = [];
for i=1:(length(pos1sol_prov)-1)
    path = pathfinder_poly([pos1sol_prov(i) pos2sol_prov(i)], [pos1sol_prov(i+1) pos2sol_prov(i+1)], P, L);
    %path = [pos1sol_prov(i) pos2sol_prov(i); pos1sol_prov(i+1) pos2sol_prov(i+1)];
    waypoints = [waypoints;path(1:(end-1),:)]; %#ok<AGROW>
end
waypoints = [waypoints;[pos1sol_prov(end) pos2sol_prov(end)]];

pos1sol = waypoints(:,1);
pos2sol = waypoints(:,2);

%-------------------------------- timestamps ----------------------------------------

t = zeros(size(pos1sol));
for i=2:length(t)
    t(i)=t(i-1)+norm([pos1sol(i)-pos1sol(i-1);pos2sol(i)-pos2sol(i-1)])/speed;
end
end