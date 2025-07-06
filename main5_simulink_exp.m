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
Ld = 10;
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

%% Compute UAR, ALOP, and Surface Coverage for each algorithm
% Each call generates a figure (saved as a PDF with the given name)

[ALOP_mst, UAR_mst, Surf_mst] = plot_and_compute_UAR(Hsize, Lsize, MagFactor, d, [pos1sol_mst'; pos2sol_mst'], pgons, npol, t_mst, speed, sprintf('Figures/UAR_MST_sim'));
