clc
clear all
close all

addpath('util')

%%
N = 356;
Hsize = 90;
Lsize = 120;
d = 25;
sigma = 20;
nsom = 12;
npol = 3;
speed = 0.4;

[P, D, L] = gen_polygons(Lsize, Hsize, sigma, npol, nsom, N);

% Obter polígonos
[pgons, npol] = recoverPolygonsFromPL(P, L);

% Planeamento
tic; [pos1sol_hex, pos2sol_hex, t_hex] = MST_plan_hex2(10,Lsize,Hsize,d*0.9,size(P, 1),npol,P,D,L,pgons,speed); comp_time_MST = toc;
tic; [pos1sol_mst, pos2sol_mst, t_mst] = MST_plan2(10,Lsize,Hsize,d*0.9,size(P, 1),npol,P,D,L,pgons,speed); comp_time_MST_hex = toc;
tic; [pos1sol, pos2sol, t] = TSP_plan(Lsize, Hsize, d, P, D, L, speed); comp_time_TSP = toc;

% OCP
rmax = 5;
pos0 = [0; 20];
yaw0 = 0;
T = 1500;
Nsteps = 80;

tic;
[pos1sol_ocp, pos2sol_ocp, t_ocp] = OCP_plan(Lsize, Hsize, 0.6*d, P, D, L, speed, pos1sol, pos2sol, max(t)*1.5);
comp_time_OCP = toc; pos1sol_ocp=pos1sol_ocp'; pos2sol_ocp=pos2sol_ocp';
%%
% Calcular UAR com a nova função
MagFactor = 4;
[ALOP_tsp, UAR_tsp, Surface_tsp] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol'; pos2sol'], pgons, npol, t, 'Figures/UAR_TSP');
[ALOP_hex, UAR_hex, Surface_hex] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_hex'; pos2sol_hex'], pgons, npol, t_hex, 'Figures/UAR_HEX');
[ALOP_mst, UAR_mst, Surface_mst] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_mst'; pos2sol_mst'], pgons, npol, t_mst, 'Figures/UAR_MST');
[ALOP_ocp, UAR_ocp, Surface_ocp] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_ocp'; pos2sol_ocp'], pgons, npol, t_ocp, 'Figures/UAR_OCP');

% Plot comparativo das trajetórias
figure
hold on
title('Trajectory')
plot(pos1sol, pos2sol, 'k-')
plot(pos1sol_hex, pos2sol_hex, 'b-')
plot(pos1sol_mst, pos2sol_mst, 'r-')
plot(pos1sol_ocp, pos2sol_ocp, 'g-')
for i = 1:max(L)
    plot(polyshape(P(L == i, 1), P(L == i, 2)));
end
print -dpdf 'Figures/overlap'

% Métricas
algorithms = {'TSP', 'MST-Hex', 'MST-square', 'OCP'};
proc_time = [comp_time_TSP, comp_time_MST_hex, comp_time_MST, comp_time_OCP];

% Tempo de processamento
figure;
bar(proc_time);
set(gca, 'YScale', 'log');
xticklabels(algorithms);
xlabel('Algorithms');
ylabel('Processing Time in s');
title('Algorithm Processing Time Comparison');
grid on;
set(gca, 'YMinorGrid', 'on', 'YMinorTick', 'on');
print -dpdf 'Figures/ProcessingTime'

% Área não coberta
uncovered_area = [Surface_tsp*UAR_tsp/100, Surface_hex*UAR_hex/100, Surface_mst*UAR_mst/100, Surface_ocp*UAR_ocp/100];
figure;
bar(uncovered_area);
xticklabels(algorithms);
xlabel('Algorithms');
ylabel('Uncovered Area in $m^2$','Interpreter','Latex');
title('Algorithm Uncovered Area Comparison');
grid on;
set(gca, 'YMinorGrid', 'on', 'YMinorTick', 'on');
print -dpdf 'Figures/UncovArea'

% Distância total percorrida
distance = [speed*t(end), speed*t_hex(end), speed*t_mst(end), speed*t_ocp(end)];
figure;
bar(distance);
xticklabels(algorithms);
xlabel('Algorithms');
ylabel('Distance in $m$','Interpreter','Latex');
title('Total path length for each planning algorithm');
grid on;
set(gca, 'YMinorGrid', 'on', 'YMinorTick', 'on');
print -dpdf 'Figures/PathLength'
