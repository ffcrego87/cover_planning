clc
clear all
close all

%%
N = 497;
% Size of the map
Hsize=90;
Lsize=120;
d = 25;
sigma = 20;% Rayon moyen des obstacles
nsom = 12;% Nombre maximum de sommets
nk = 0;
npol = 2;
speed = 0.4;

[P,D,L]=gen_polygons(Lsize,Hsize,sigma,npol,nsom,N);

[pos1sol_mst,pos2sol_mst,t_mst] = MST_plan(Lsize,Hsize,d,P,D,L,speed);
[pos1sol_hex,pos2sol_hex,t_hex] = MST_plan_hex(Lsize,Hsize,d,P,D,L,speed);
[pos1sol,pos2sol,t] = TSP_plan(Lsize,Hsize,d,P,D,L,speed);

% Problem settings OCP
rmax = 5;
pos0 = [0;20];
yaw0 = 0;
T = 1500;  % final time

% OCP settings
N = 80; % number of control intervals

[pos1sol_ocp,pos2sol_ocp,t_ocp] = OCP_plan(Lsize,Hsize,d,P,D,L,speed,T,rmax,pos0,yaw0,N);

%%
%Plotting
figure
hold on
plot(t,pos1sol,'b-')
plot(t,pos2sol,'r-')
xlabel('Time (s)')

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

%%

% disc=1;
% 
% xrange = 0:disc:Lsize;
% yrange = 0:disc:Hsize;
% 
% Cellpos1 = kron(xrange,ones(size(yrange)));
% Cellpos2 = kron(ones(size(xrange)),yrange);
% 
% idxs_to_remove = [];
% 
% 
% figure
% hold on
% 
% length(Cellpos1)
% 
% 
% plot(pos1sol_ocp,pos2sol_ocp,'k-')
% 
% for i=1:max(L)
%     Pi1 = P(L==i,:);
%     Pi = [Pi1;Pi1(1,:)];
%     plot(polyshape(P(L==i,1),P(L==i,2)));
%     for k=1:length(Cellpos1)
%         if inpolygon(Cellpos1(k),Cellpos2(k),Pi(:,1),Pi(:,2))
%             idxs_to_remove = [idxs_to_remove;k];
%         end
%     end
% end
% 
% Cellpos1(idxs_to_remove) = [];
% Cellpos2(idxs_to_remove) = [];
% 
% idxs_to_remove = [];
% 
% for k=1:length(Cellpos1)
%     for l = 1:length(pos1sol_ocp)
%         if norm([Cellpos1(k);Cellpos2(k)]-[pos1sol_ocp(l);pos2sol_ocp(l)])<(d/2)
%             idxs_to_remove = [idxs_to_remove;k];
%             continue;
%         end
%     end
% end
% 
% Cellpos1(idxs_to_remove) = [];
% Cellpos2(idxs_to_remove) = [];
% 
% uncovered_area = length(Cellpos1)*disc^2;
% 
% plot(Cellpos1,Cellpos2,'k.')




