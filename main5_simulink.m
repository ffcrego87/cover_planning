clc;
clear all;
close all;

addpath('util');

%% Fixed Parameters
nSamples = 10;  % Number of samples (editable)
N = 4;
Hsize = 90;
Lsize = 120;
d = 25;
sigma = 20;
nsom = 12;
npol_fixed = 3;  % Number of polygons
speed = 0.4;
MagFactor = 4;

% sim parameters
K = 0.5;
V = 2.5;
T = 5;
KP = 1.2;
KD = 3;
Ld = 20;
scale = 10;

% Create folder for saving figures, if it doesn't exist
if ~exist('Figures', 'dir')
    mkdir('Figures');
end

% Matrices to store results [sample x algorithm]
% Algorithm order: 1 - TSP, 2 - MST-Hex, 3 - MST-Square, 4 - OCP, 5 -
% Boustrophedon

fprintf('Running sims... \n');

% Generate polygons
[P, D, L] = gen_polygons(Lsize, Hsize, sigma, npol_fixed, nsom, N);

% Recover polygons
[pgons, npol] = recoverPolygonsFromPL(P, L);

%% Trajectory planning for each algorithm

% --- TSP ---
tic;
[pos1sol, pos2sol, t] = TSP_plan(Lsize, Hsize, d, P, D, L, speed);
comp_time_TSP = toc;

pos1sol_prev = pos1sol;
pos2sol_prev = pos2sol;

waypoints = scale*[pos1sol, pos2sol];
x0 = waypoints(1,1);
y0 = waypoints(1,2);
psi_0=0;
out=sim('planning.slx');
pos1sol=out.pos(:,1)/scale;
pos2sol=out.pos(:,2)/scale;

% --- MST - Hexagonal ---
tic;
[pos1sol_hex, pos2sol_hex, t_hex] = MST_plan_hex2(10, Lsize, Hsize, d*0.9, size(P,1), npol, P, D, L, pgons, speed);
comp_time_MST_hex = toc;

waypoints = scale*[pos1sol_hex, pos2sol_hex];
x0 = waypoints(1,1);
y0 = waypoints(1,2);
psi_0=0;
out=sim('planning.slx');
pos1sol_hex=out.pos(:,1)/scale;
pos2sol_hex=out.pos(:,2)/scale;

% --- MST - Square ---
tic;
[pos1sol_mst, pos2sol_mst, t_mst] = MST_plan2(10, Lsize, Hsize, d*0.9, size(P,1), npol, P, D, L, pgons, speed);
comp_time_MST = toc;

waypoints = scale*[pos1sol_mst, pos2sol_mst];
x0 = waypoints(1,1);
y0 = waypoints(1,2);
psi_0=0;
out=sim('planning.slx');
pos1sol_mst=out.pos(:,1)/scale;
pos2sol_mst=out.pos(:,2)/scale;

% --- MST - Square ---
tic;
[pos1sol_bous, pos2sol_bous, t_bous] = Boustro_backtrack_plan(Lsize, Hsize, d, P, D, L, speed);
comp_time_bous = toc;

waypoints = scale*[pos1sol_bous, pos2sol_bous];
x0 = waypoints(1,1);
y0 = waypoints(1,2);
psi_0=0;
out=sim('planning.slx');
pos1sol_bous=out.pos(:,1)/scale;
pos2sol_bous=out.pos(:,2)/scale;

%%
% --- OCP ---
tic;
[pos1sol_ocp, pos2sol_ocp, t_ocp] = OCP_plan(Lsize, Hsize, 0.55*d, P, D, L, speed, pos1sol_prev, pos2sol_prev, max(t)*1.5);
comp_time_OCP = toc;
pos1sol_ocp = pos1sol_ocp';
pos2sol_ocp = pos2sol_ocp';

waypoints = scale*[pos1sol_ocp pos2sol_ocp];
x0 = waypoints(1,1);
y0 = waypoints(1,2);
psi_0=0;
out=sim('planning.slx');
pos1sol_ocp=out.pos(:,1)/scale;
pos2sol_ocp=out.pos(:,2)/scale;

%% Compute UAR, ALOP, and Surface Coverage for each algorithm
% Each call generates a figure (saved as a PDF with the given name)

[ALOP_tsp, UAR_tsp, Surf_tsp] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol'; pos2sol'], pgons, npol, t, speed, sprintf('Figures/UAR_TSP_sim'));
[ALOP_hex, UAR_hex, Surf_hex] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_hex'; pos2sol_hex'], pgons, npol, t_hex, speed, sprintf('Figures/UAR_HEX_sim'));
[ALOP_mst, UAR_mst, Surf_mst] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_mst'; pos2sol_mst'], pgons, npol, t_mst, speed, sprintf('Figures/UAR_MST_sim'));
[ALOP_ocp, UAR_ocp, Surf_ocp] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_ocp'; pos2sol_ocp'], pgons, npol, t_ocp, speed, sprintf('Figures/UAR_OCP_sim'));
[ALOP_bous, UAR_bous, Surf_bous] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_bous'; pos2sol_bous'], pgons, npol, t_bous, speed, sprintf('Figures/UAR_bous_sim'));

% Store metrics for each algorithm
ALOP_mat     = [ALOP_tsp,   ALOP_hex,   ALOP_mst,   ALOP_ocp,   ALOP_bous];
UAR_mat      = [UAR_tsp,    UAR_hex,    UAR_mst,    UAR_ocp,    UAR_bous];
compTime_mat = [comp_time_TSP, comp_time_MST_hex, comp_time_MST, comp_time_OCP, comp_time_bous];
surf_mat     = [Surf_tsp, Surf_hex, Surf_mst, Surf_ocp, Surf_bous];

algorithms = {'TSP', 'MST-Hex', 'MST-Square', 'OCP', 'Boustrophedon'};

%% Bar graphs with mean and standard deviation

% ALOP Plot
figure;
bar(ALOP_mat);
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithms');
ylabel('ALOP');
title('ALOP');
grid on;
print -dpdf 'Figures/ALOP_Sim';

% Computation Time Plot (with logarithmic scale)
figure;
bar(compTime_mat);
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithms');
ylabel('Computation Time (s)');
title('Average Computation Time');
grid on;
set(gca, 'YScale', 'log');  % Define y-axis in logarithmic scale
print -dpdf 'Figures/ComputationTime_Sim';

% UAR Plot
figure;
bar(UAR_mat);
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithms');
ylabel('UAR (%)');
title('Average UAR');
grid on;
print -dpdf 'Figures/UAR_Sim';

% Surface Coverage Plot
figure;
bar(surf_mat);
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithms');
ylabel('Surface Coverage (%)');
title('Average Surface Coverage with Standard Deviation');
grid on;
print -dpdf 'Figures/SurfaceCoverage_Sim';

%% Build summary table and save to file  -------------------------------
T = table( algorithms(:),                                   ... % 1st col
           ALOP_mat(:),                     ... % ALOP
           compTime_mat(:),              ... % comp-time
           UAR_mat(:),                      ... % UAR
           surf_mat(:),                    ... % coverage
           'VariableNames', ...
          {'Algorithm' ,  'ALOP', ...
           'CompTime', ...
           'UAR', ...
           'Surf'} );

% choose your preferred format
outDir = 'Results';            % create folder if it doesn’t exist
if ~exist(outDir,'dir'), mkdir(outDir); end

writetable(T, fullfile(outDir,'metrics_summary_sim.csv'));   % CSV
% writetable(T, fullfile(outDir,'metrics_summary.xlsx')); % Excel alternative