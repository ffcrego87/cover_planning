clc;
clear all;
close all;

addpath('util');

%% Fixed Parameters
nSamples = 10;  % Number of samples (editable)
N = 356;
Hsize = 90;
Lsize = 120;
d = 25;
sigma = 20;
nsom = 12;
npol_fixed = 3;  % Number of polygons
speed = 0.4;
MagFactor = 4;

% Create folder for saving figures, if it doesn't exist
if ~exist('Figures', 'dir')
    mkdir('Figures');
end

% Matrices to store results [sample x algorithm]
% Algorithm order: 1 - TSP, 2 - MST-Hex, 3 - MST-Square, 4 - OCP
ALOP_mat     = zeros(nSamples, 4);
UAR_mat      = zeros(nSamples, 4);
compTime_mat = zeros(nSamples, 4);
surf_mat     = zeros(nSamples, 4);  % New: Surface coverage metric

%% Sampling Loop
for sample = []
    fprintf('Running sample %d...\n', sample);
    
    % Generate polygons
    [P, D, L] = gen_polygons(Lsize, Hsize, sigma, npol_fixed, nsom, N);
    
    % Recover polygons
    [pgons, npol] = recoverPolygonsFromPL(P, L);
    
    %% Trajectory planning for each algorithm
    
    % --- TSP ---
    tic;
    [pos1sol, pos2sol, t] = TSP_plan(Lsize, Hsize, d, P, D, L, speed);
    comp_time_TSP = toc;
    
    % --- MST - Hexagonal ---
    tic;
    [pos1sol_hex, pos2sol_hex, t_hex] = MST_plan_hex2(10, Lsize, Hsize, d*0.9, size(P,1), npol, P, D, L, pgons, speed);
    comp_time_MST_hex = toc;
    
    % --- MST - Square ---
    tic;
    [pos1sol_mst, pos2sol_mst, t_mst] = MST_plan2(10, Lsize, Hsize, d*0.9, size(P,1), npol, P, D, L, pgons, speed);
    comp_time_MST = toc;
    
    % --- OCP ---
    rmax = 5;
    pos0 = [0; 20];
    yaw0 = 0;
    T_val = 1500;
    Nsteps = 80;
    tic;
    [pos1sol_ocp, pos2sol_ocp, t_ocp] = OCP_plan(Lsize, Hsize, 0.6*d, P, D, L, speed, pos1sol, pos2sol, max(t)*1.5);
    comp_time_OCP = toc;
    pos1sol_ocp = pos1sol_ocp';
    pos2sol_ocp = pos2sol_ocp';
    
    %% Compute UAR, ALOP, and Surface Coverage for each algorithm
    % Each call generates a figure (saved as a PDF with the given name)
    
    [ALOP_tsp, UAR_tsp, Surf_tsp] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol'; pos2sol'], pgons, npol, t, speed, sprintf('Figures/UAR_TSP_sample%d', sample));
    [ALOP_hex, UAR_hex, Surf_hex] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_hex'; pos2sol_hex'], pgons, npol, t_hex, speed, sprintf('Figures/UAR_HEX_sample%d', sample));
    [ALOP_mst, UAR_mst, Surf_mst] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_mst'; pos2sol_mst'], pgons, npol, t_mst, speed, sprintf('Figures/UAR_MST_sample%d', sample));
    [ALOP_ocp, UAR_ocp, Surf_ocp] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_ocp'; pos2sol_ocp'], pgons, npol, t_ocp, speed, sprintf('Figures/UAR_OCP_sample%d', sample));
    
    % Store metrics for each algorithm
    ALOP_mat(sample,:)     = [ALOP_tsp,   ALOP_hex,   ALOP_mst,   ALOP_ocp];
    UAR_mat(sample,:)      = [UAR_tsp,    UAR_hex,    UAR_mst,    UAR_ocp];
    compTime_mat(sample,:) = [comp_time_TSP, comp_time_MST_hex, comp_time_MST, comp_time_OCP];
    surf_mat(sample,:)     = [Surf_tsp, Surf_hex, Surf_mst, Surf_ocp];  % New metric
end

%% Compute means and standard deviations
avg_ALOP     = mean(ALOP_mat);
std_ALOP     = std(ALOP_mat);

avg_compTime = mean(compTime_mat);
std_compTime = std(compTime_mat);

avg_UAR      = mean(UAR_mat);
std_UAR      = std(UAR_mat);

avg_surf     = mean(surf_mat);   % New
std_surf     = std(surf_mat);    % New

algorithms = {'TSP', 'MST-Hex', 'MST-Square', 'OCP'};

%% Bar graphs with mean and standard deviation

% ALOP Plot
figure;
bar(avg_ALOP);
hold on;
errorbar(1:length(avg_ALOP), avg_ALOP, std_ALOP, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithms');
ylabel('ALOP');
title('Average ALOP with Standard Deviation');
grid on;
print -dpdf 'Figures/ALOP_Comparison';

% Computation Time Plot (with logarithmic scale)
figure;
bar(avg_compTime);
hold on;
errorbar(1:length(avg_compTime), avg_compTime, std_compTime, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithms');
ylabel('Computation Time (s)');
title('Average Computation Time with Standard Deviation');
grid on;
set(gca, 'YScale', 'log');  % Define y-axis in logarithmic scale
print -dpdf 'Figures/ComputationTime_Comparison';

% UAR Plot
figure;
bar(avg_UAR);
hold on;
errorbar(1:length(avg_UAR), avg_UAR, std_UAR, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithms');
ylabel('UAR (%)');
title('Average UAR with Standard Deviation');
grid on;
print -dpdf 'Figures/UAR_Comparison';

% Surface Coverage Plot (new)
figure;
bar(avg_surf);
hold on;
errorbar(1:length(avg_surf), avg_surf, std_surf, 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
set(gca, 'XTickLabel', algorithms);
xlabel('Algorithms');
ylabel('Surface Coverage (%)');
title('Average Surface Coverage with Standard Deviation');
grid on;
print -dpdf 'Figures/SurfaceCoverage_Comparison';
