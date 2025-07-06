function [pos1sol,pos2sol,t] = Boustro_plan(Lsize,Hsize,d,P,~,L,speed)
% BOUSTRO_PLAN  Very simple boustrophedon path generator on a regular grid.
%
%   [xPath,yPath,t] = Boustro_plan(Lsize,Hsize,d,P,~,L,speed)
%   - Lsize,Hsize   : size of the bounding box (m)
%   - d             : tool width (m)              (we sample at d/2)
%   - P , L         : vertices & labels of polygonal obstacles
%   - speed         : robot linear speed (m/s)
%
%   OUTPUT:
%     pos1sol  = column vector of x-coordinates   [m]
%     pos2sol  = column vector of y-coordinates   [m]
%     t        = time stamps assuming constant speed [s]

% ------------------------------------------------------------
% 1. Build a half-resolution grid and remove points inside obstacles
% ------------------------------------------------------------
step = d/2;                                      % lateral sample
[X,Y] = meshgrid(0:step:Lsize , 0:step:Hsize);
v     = [X(:) Y(:)];

obstacles = cell(max(L),1);
for i = 1:max(L), obstacles{i} = P(L==i,:); end

for i = 1:numel(obstacles)
    obs = obstacles{i};
    [in,on] = inpolygon(v(:,1),v(:,2), obs(:,1),obs(:,2));
    v = v(~(in|on),:);                           % keep only free cells
end

% ------------------------------------------------------------
% 2. Order points in a boustrophedon sequence
% ------------------------------------------------------------
% 2.a  group by unique Y rows (tolerance = step/10 to absorb fp error)
[rowsY,~,rowIdx] = uniquetol(v(:,2), step/100);   % ascending Y

path = zeros(size(v)); k = 0;
for r = 1:numel(rowsY)
    ptsRow  = v(rowIdx==r,:);                    % all pts on this Y
    if mod(r,2)==1                               % odd row  -> left → right
        disp('ascend')
        ptsRow = sortrows(ptsRow,1,'ascend');
    else                                         % even row -> right → left
        disp('descend')
        ptsRow = sortrows(ptsRow,1,'descend');
    end
    path(k+(1:size(ptsRow,1)),:) = ptsRow;       % append to master path
    k = k + size(ptsRow,1);
end
path = path(1:k,:);                              % trim unused prealloc

pos1sol = path(:,1);                             % x-coords
pos2sol = path(:,2);                             % y-coords

% ------------------------------------------------------------
% 3. Build time vector (constant-speed travel between way-points)
% ------------------------------------------------------------
t = zeros(size(pos1sol));
for i = 2:numel(t)
    segLen   = hypot(pos1sol(i)-pos1sol(i-1), pos2sol(i)-pos2sol(i-1));
    t(i)     = t(i-1) + segLen/speed;
end
end