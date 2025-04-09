%% Clear workspace and add utility functions path
clc;
clear all;
close all;
addpath('util');

%% Map and planning parameters
Lsize = 600;        % Map width (units)
Hsize = 450;        % Map height (units)
d = 20;             % Safety distance (or resolution)
speed = 0.4;        % Speed parameter

%% List of images (assumed to be in the "Maps" folder)
images = {'snazzy-image.png', 'snazzy-image2.png', 'snazzy-image3.png'};
num_maps = length(images);

%% Initialize an array to store results (total distance traveled)
dist_mst = zeros(num_maps, 1);

%% Loop through each image for processing
for m = 1:num_maps
    fprintf('Processing image %d of %d...\n', m, num_maps);
    image_file = fullfile('Maps', images{m});
    
    % Read the original image
    RGB = imread(image_file);
    
    % Create a mask using a helper function
    [BW, maskedRGB] = createMask(RGB);
    
    % Extract the contour of the region using contour_path
    % Syntax assumed: [xContour, yContour] = contour_path(dim, value, BW)
    [xContour, yContour] = contour_path(size(BW,1), 1, BW);
    
    % Build the polygon from the contour and close the polygon
    % so that the first point equals the last point.
    polygon_out = [xContour' xContour(1); yContour' yContour(1)]';
    
    % Simplify the polygon if necessary
    polygon_out = simplify_pgon(polygon_out);
    
    %% Trajectory planning using MST_plan_out
    % Call the planning function with the specified parameters.
    % The empty arrays ([]) indicate unused parameters.
    [xSol, ySol, timeSol] = MST_plan_out(Lsize, Hsize, d, [], [], [], polygon_out, speed);
    
    % Calculate the total distance traveled along the trajectory
    total_distance = sum(sqrt(diff(xSol).^2 + diff(ySol).^2));
    dist_mst(m) = total_distance;
    
    %% Plot the original image as background, region polygon, and trajectory
    figure;
    
    % Display the original image as background
    imshow(RGB);
    hold on;
    
    % Optionally, adjust the axis to match the image size if necessary
    % axis([1 size(RGB,2) 1 size(RGB,1)]);
    
    % Plot the region polygon with transparency so that the background shows through
    regionShape = polyshape(polygon_out(:,1), size(RGB,1)-polygon_out(:,2), 'Simplify', false);
    plot(regionShape, 'FaceColor', [0.9 0.9 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5);
    
    % Plot the MST trajectory on top of the image
    plot(xSol, size(RGB,1)-ySol, 'b-', 'LineWidth', 2);
    
    title(sprintf('Map %d - MST Trajectory', m));
    legend('Region', 'MST Trajectory', 'Location', 'best');
    grid on;
    hold off;
    
    % Option: Save the figure as a PDF file (uncomment the next line if desired)
    print('-dpdf', fullfile('Figures', sprintf('Map_%d_MST_Trajectory.pdf', m)));
end

%% Display total distances for each map (optional)
fprintf('\nTotal distances for MST OUT trajectory:\n');
for m = 1:num_maps
    fprintf('Map %d: %.2f units\n', m, dist_mst(m));
end
