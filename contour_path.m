function [pos1sol,pos2sol] = contour_path(Hsize,scale,BW)

Img = 255*uint8(BW);

%%
BW3 = bwmorph(Img,'remove');
BW3 = bwskel(BW3,'MinBranchLength',20);% Almost unnecessary except for some few cases, such as Map N°252 or 221 (requiring MinBranchLength>=8)

%% Find parametrized path
Mpath = BW3;
no=10;
ymax = scale*Hsize;
npath = sum(sum(Mpath));
path = zeros(npath,2);
[yp,xp] = find(Mpath,1);
path(1,:) = [xp,ymax-yp]/scale;
Mpath(yp,xp) = false;
for j=2:npath
    [dy,dx] = find(Mpath(max(yp-no,1):min(yp+no,size(Mpath,1)),max(xp-no,1):min(xp+no,size(Mpath,2))));
    if isempty(dx)
        if npath-j>0 % strange
            warning('Forgot %d points, probably a matter of connexity of overlapping\nTry increase MinBranchLength in bwskel\n',npath-j)
        end
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
pos1sol=path(1:(j-1),1);
pos2sol=path(1:(j-1),2);
end