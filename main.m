clc
clear all
close all

addpath('util')

%%
N = 356;%352;356
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


tic; [pos1sol_mst,pos2sol_mst,t_mst] = MST_plan(Lsize,Hsize,d,P,D,L,speed); comp_time_MST=toc;
tic; [pos1sol_hex,pos2sol_hex,t_hex] = MST_plan_hex(Lsize,Hsize,d,P,D,L,speed); comp_time_MST_hex=toc;
tic; [pos1sol,pos2sol,t] = TSP_plan(Lsize,Hsize,d,P,D,L,speed); comp_time_TSP=toc;

% Problem settings OCP
rmax = 5;
pos0 = [0;20];
yaw0 = 0;
T = 1500;  % final time

% OCP settings
N = 80; % number of control intervals
%%
tic; [pos1sol_ocp,pos2sol_ocp,t_ocp] = OCP_plan(Lsize,Hsize,0.8*d,P,D,L,speed,pos1sol,pos2sol,max(t)*1.5); comp_time_OCP=toc;

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
print -dpdf 'Figures/overlap'

uncovered_area_tsp = uncovered_area_comp(pos1sol,pos2sol,d,Lsize,Hsize,P,D,L);title('Trajectory with TSP');
print -dpdf 'Figures/TrajectoryTSP'
uncovered_area_mst = uncovered_area_comp(pos1sol_mst,pos2sol_mst,d,Lsize,Hsize,P,D,L);title('Trajectory with MST Square Configuration');
print -dpdf 'Figures/TrajectoryMST'
uncovered_area_hex = uncovered_area_comp(pos1sol_hex,pos2sol_hex,d,Lsize,Hsize,P,D,L);title('Trajectory with MST Hexagonal Configuration');
print -dpdf 'Figures/TrajectoryHEX'
uncovered_area_ocp = uncovered_area_comp(pos1sol_ocp,pos2sol_ocp,d,Lsize,Hsize,P,D,L);title('Trajectory with OCP');
print -dpdf 'Figures/TrajectoryOCP'


% Data
algorithms = {'TSP', 'MST-Hex', 'MST-square', 'OCP'};
distance = [comp_time_TSP, comp_time_MST_hex, comp_time_MST, comp_time_OCP];

% Create bar plot
figure;
bar(distance);
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
xticklabels(algorithms); % Set x-axis labels
xlabel('Algorithms');
ylabel('Processing Time in s');
title('Algorithm Processing Time Comparison');

% Improve grid and appearance
grid on;
set(gca, 'YMinorGrid', 'on', 'YMinorTick', 'on');
print -dpdf 'Figures/ProcessingTime'



% Data
uncovered_area = [uncovered_area_tsp, uncovered_area_hex, uncovered_area_mst, uncovered_area_ocp];

% Create bar plot
figure;
bar(uncovered_area);
xticklabels(algorithms); % Set x-axis labels
xlabel('Algorithms');
ylabel('Uncovered Area in $m^2$','Interpreter','Latex');
title('Algorithm Uncovered Area Comparison');

% Improve grid and appearance
grid on;
set(gca, 'YMinorGrid', 'on', 'YMinorTick', 'on');
print -dpdf 'Figures/UncovArea'



% Data
distance = [speed*t(end), speed*t_hex(end), speed*t_mst(end), speed*t_ocp(end)];

% Create bar plot
figure;
bar(distance);
xticklabels(algorithms); % Set x-axis labels
xlabel('Algorithms');
ylabel('Distance in $m$','Interpreter','Latex');
title('Total path length for each planning algorithm');

% Improve grid and appearance
grid on;
set(gca, 'YMinorGrid', 'on', 'YMinorTick', 'on');
print -dpdf 'Figures/PathLength'




