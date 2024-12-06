clc
clear all
close all

%%
N = 350;
% Size of the map
Hsize=450;
Lsize=600;
d = 20;
sigma = 100;% Rayon moyen des obstacles
nsom = 12;% Nombre maximum de sommets
nk = 0;
npol = 1;
speed = 0.4;

RGB = imread('snazzy-image3.png');
[BW,maskedRGBImage] = createMask(RGB);
[pos1sol,pos2sol] = contour_path(size(BW,1),1,BW);
Pout = [pos1sol' pos1sol(1); pos2sol' pos2sol(1)]';

Pout = simplify_pgon(Pout);

%%

[pos1sol_mst,pos2sol_mst,t_mst] = MST_plan_out(Lsize,Hsize,d,[],[],[],Pout,speed);
[pos1sol_hex,pos2sol_hex,t_hex] = MST_plan_hex_out(Lsize,Hsize,d,[],[],[],Pout,speed);
%[pos1sol_tsp,pos2sol_tsp,t_tsp] = TSP_plan2_out(Lsize,Hsize,d,[],[],[],Pout,speed);
%%
figure
hold on
title('Trajectory')
% Remove duplicate vertices
plot(polyshape(Pout(:,1),Pout(:,2),'Simplify',false))
plot(pos1sol_mst,pos2sol_mst,'b-')

figure
hold on
title('Trajectory')
% Remove duplicate vertices
plot(polyshape(Pout(:,1),Pout(:,2),'Simplify',false))
plot(pos1sol_hex,pos2sol_hex,'r-')

% figure
% hold on
% title('Trajectory')
% % Remove duplicate vertices
% plot(polyshape(Pout(:,1),Pout(:,2),'Simplify',false))
% plot(pos1sol_tsp,pos2sol_tsp,'k-')





