RGB = imread('snazzy-image3.png');
[BW,maskedRGBImage] = createMask(RGB);
[pos1sol,pos2sol] = contour_path(size(BW,1),1,BW);
pgon = polyshape([pos1sol' pos1sol(1)],[pos2sol' pos2sol(1)]);
plot(pgon)