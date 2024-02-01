clc
clear all
close all

%%
N = 350;
% Size of the map
Hsize=90;
Lsize=120;
d = 25;
sigma = 20;% Rayon moyen des obstacles
nsom = 12;% Nombre maximum de sommets
nk = 0;
npol = 3;
speed = 0.4;

[P,D,L]=gen_polygons(Lsize,Hsize,sigma,npol,nsom,N);

[pos1sol_mst,pos2sol_mst,t_mst] = MST_plan(Lsize,Hsize,d,P,D,L,speed);
[pos1sol_hex,pos2sol_hex,t_hex] = MST_plan_hex(Lsize,Hsize,d,P,D,L,speed);
[pos1sol,pos2sol,t] = TSP_plan2(Lsize,Hsize,d,P,D,L,speed);

% Problem settings OCP
rmax = 5;
pos0 = [0;20];
yaw0 = 0;
T = 1500;  % final time

% OCP settings
N = 80; % number of control intervals
%%
[pos1sol_ocp,pos2sol_ocp,t_ocp] = OCP_plan4(Lsize,Hsize,0.8*d,P,D,L,speed,pos1sol,pos2sol,max(t)*1.5);

%%
%Plotting

figure
hold on
title('Trajectory')
plot(pos1sol,pos2sol,'k-')
plot(pos1sol_hex,pos2sol_hex,'b-')
plot(pos1sol_mst,pos2sol_mst,'r-')
plot(pos1sol_ocp,pos2sol_ocp,'g-')
for i=1:max(L)
    plot(polyshape(P(L==i,1),P(L==i,2)));
end

uncovered_area_tsp = uncovered_area_comp(pos1sol,pos2sol,d,Lsize,Hsize,P,D,L);
uncovered_area_mst = uncovered_area_comp(pos1sol_mst,pos2sol_mst,d,Lsize,Hsize,P,D,L);
uncovered_area_hex = uncovered_area_comp(pos1sol_hex,pos2sol_hex,d,Lsize,Hsize,P,D,L);
uncovered_area_ocp = uncovered_area_comp(pos1sol_ocp,pos2sol_ocp,d,Lsize,Hsize,P,D,L);
% uncovered_area = length(Cellpos1)*disc^2;
% 
% plot(Cellpos1,Cellpos2,'k.')




