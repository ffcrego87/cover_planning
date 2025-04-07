images = {'snazzy-image.png', 'snazzy-image2.png', 'snazzy-image3.png'};
num_maps = length(images);
% Initialize arrays for results, e.g., for MST
dist_mst = zeros(num_maps, 1);
covered_mst = zeros(num_maps, 1);
% Repeat for other methods (hexagonal MST, TSP, OCP)

for m = 1:num_maps
    image_file = images{m};
    P = detect_polygons_from_image(image_file, Lsize, Hsize);
    % Generate points in free space
    num_candidates = 1000;
    x_cand = rand(num_candidates, 1) * Lsize;
    y_cand = rand(num_candidates, 1) * Hsize;
    inside = false(num_candidates, 1);
    for k = 1:length(P)
        inside = inside | inpolygon(x_cand, y_cand, P{k}(:,1), P{k}(:,2));
    end
    free_points = [x_cand(~inside), y_cand(~inside)];
    if size(free_points, 1) < N
        error('Not enough free space points for N=%d', N);
    end
    idx = randperm(size(free_points, 1), N);
    points = free_points(idx, :);
    % Run path planning, e.g., for MST
    [pos1sol_mst, pos2sol_mst, t_mst] = MST_plan(points, P, d, ...); % Adjust parameters
    dist_mst(m) = sum(sqrt(diff(pos1sol_mst).^2 + diff(pos2sol_mst).^2));
    % Compute covered area
    uncovered = uncovered_area_comp(pos1sol_mst, pos2sol_mst, t_mst, P, d, Lsize, Hsize);
    obstacle_area = sum(cellfun(@(p) polyarea(p(:,1), p(:,2)), P));
    total_area = Lsize * Hsize;
    free_space_area = total_area - obstacle_area;
    covered_mst(m) = free_space_area - uncovered;
    % Repeat for other methods
end

% Compute averages
avg_dist_mst = mean(dist_mst);
avg_covered_mst = mean(covered_mst);
% Display results in a table for clarity
results = table({'MST'; 'Hexagonal MST'; 'TSP'; 'OCP'}, ...
    [avg_dist_mst; avg_dist_hex; avg_dist_tsp; avg_dist_ocp], ...
    [avg_covered_mst; avg_covered_hex; avg_covered_tsp; avg_covered_ocp], ...
    'VariableNames', {'Method', 'AvgDistance', 'AvgCoveredArea'});
disp(results);