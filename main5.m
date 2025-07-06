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
% Algorithm order: 1 - TSP, 2 - MST-Hex, 3 - MST-Square, 4 - OCP, 5 -
% Boustrophedon
ALOP_mat     = zeros(nSamples, 5);
UAR_mat      = zeros(nSamples, 5);
compTime_mat = zeros(nSamples, 5);
surf_mat     = zeros(nSamples, 5);

%% Sampling Loop
smpidx=1;
for sample = [1,2,3,4,5,6,8,9,10,11]
    fprintf('Running sample %d...\n', sample);
    
    % Generate polygons
    [P, D, L] = gen_polygons(Lsize, Hsize, sigma, npol_fixed, nsom, sample);
    
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
    
    % --- MST - Square ---
    tic;
    [pos1sol_bous, pos2sol_bous, t_bous] = Boustro_backtrack_plan(Lsize, Hsize, d, P, D, L, speed);
    comp_time_bous = toc;
    
    % --- OCP ---
    rmax = 5;
    pos0 = [0; 20];
    yaw0 = 0;
    T_val = 1500;
    Nsteps = 80;
    tic;
    [pos1sol_ocp, pos2sol_ocp, t_ocp] = OCP_plan(Lsize, Hsize, 0.55*d, P, D, L, speed, pos1sol, pos2sol, max(t)*1.5);
    comp_time_OCP = toc;
    pos1sol_ocp = pos1sol_ocp';
    pos2sol_ocp = pos2sol_ocp';
    
    %% Compute UAR, ALOP, and Surface Coverage for each algorithm
    % Each call generates a figure (saved as a PDF with the given name)
    
    [ALOP_tsp, UAR_tsp, Surf_tsp] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol'; pos2sol'], pgons, npol, t, speed, sprintf('Figures/UAR_TSP_sample%d', sample));
    [ALOP_hex, UAR_hex, Surf_hex] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_hex'; pos2sol_hex'], pgons, npol, t_hex, speed, sprintf('Figures/UAR_HEX_sample%d', sample));
    [ALOP_mst, UAR_mst, Surf_mst] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_mst'; pos2sol_mst'], pgons, npol, t_mst, speed, sprintf('Figures/UAR_MST_sample%d', sample));
    [ALOP_ocp, UAR_ocp, Surf_ocp] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_ocp'; pos2sol_ocp'], pgons, npol, t_ocp, speed, sprintf('Figures/UAR_OCP_sample%d', sample));
    [ALOP_bous, UAR_bous, Surf_bous] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_bous'; pos2sol_bous'], pgons, npol, t_bous, speed, sprintf('Figures/UAR_bous_sample%d', sample));
    
    % Store metrics for each algorithm
    ALOP_mat(smpidx,:)     = [ALOP_tsp,   ALOP_hex,   ALOP_mst,   ALOP_ocp,   ALOP_bous];
    UAR_mat(smpidx,:)      = [UAR_tsp,    UAR_hex,    UAR_mst,    UAR_ocp,    UAR_bous];
    compTime_mat(smpidx,:) = [comp_time_TSP, comp_time_MST_hex, comp_time_MST, comp_time_OCP, comp_time_bous];
    surf_mat(smpidx,:)     = [Surf_tsp, Surf_hex, Surf_mst, Surf_ocp, Surf_bous];
    smpidx=smpidx+1;
end

%% Compute means and standard deviations
avg_ALOP     = mean(ALOP_mat);
std_ALOP     = std(ALOP_mat);

avg_compTime = mean(compTime_mat);
std_compTime = std(compTime_mat);

avg_UAR      = mean(UAR_mat);
std_UAR      = std(UAR_mat);

avg_surf     = mean(surf_mat);   
std_surf     = std(surf_mat);    

algorithms = {'TSP', 'MST-Hex', 'MST-Square', 'OCP', 'Boustrophedon'};

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

%% Build summary table and save to file  -------------------------------
T = table( algorithms(:),                                   ... % 1st col
           avg_ALOP(:),   std_ALOP(:),                     ... % ALOP
           avg_compTime(:), std_compTime(:),               ... % comp-time
           avg_UAR(:),    std_UAR(:),                      ... % UAR
           avg_surf(:),   std_surf(:),                     ... % coverage
           'VariableNames', ...
          {'Algorithm' ,  'ALOP_mean','ALOP_std', ...
           'CompTime_mean','CompTime_std', ...
           'UAR_mean','UAR_std', ...
           'Surf_mean','Surf_std'} );

% choose your preferred format
outDir = 'Results';            % create folder if it doesn’t exist
if ~exist(outDir,'dir'), mkdir(outDir); end

writetable(T, fullfile(outDir,'metrics_summary.csv'));   % CSV
