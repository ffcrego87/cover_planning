function Im = roundCorners(Im0,r)
r=round(r);
Im = Im0;
Im(1,r:end-r)=0;
Im(end,r:end-r)=0;
Im(r:end-r,1)=0;
Im(r:end-r,end)=0;
for t=0:0.001:pi/2
    Im(round(1+r-r*cos(t)),round(1+r-r*sin(t)))=0;
    Im(round(1+r-r*cos(t)),round(end-1-r+r*sin(t)))=0;
    Im(round(end-1-r+r*cos(t)),round(1+r-r*sin(t)))=0;
    Im(round(end-1-r+r*cos(t)),round(end-1-r+r*sin(t)))=0;
end
end