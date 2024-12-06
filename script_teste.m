clear all;
close all;
clc;

% Create a binary image with the given shape boundary
imageSize = [14, 12];
binaryImage = zeros(imageSize);

% Define vertices of the shape in clockwise direction
shapeVertices = [
    1, 4; 1, 6; 3, 6; 3, 10; 5, 10; 5, 6; 10, 6; 10, 9;
    8, 9; 8, 11; 12, 11; 12, 9; 14, 9; 14, 6; 11, 6;
    11, 3; 14, 3; 14, 1; 10, 1; 10, 3; 7, 3; 7, 1; 4, 1;
    4, 3; 3, 3; 3, 4; 1, 4
];

% Set the shape boundary in the binary image
shapeBoundary = poly2mask(shapeVertices(:, 2), shapeVertices(:, 1), imageSize(1), imageSize(2));
binaryImage(shapeBoundary) = 1;
% Find contours in the binary image
contours = bwboundaries(binaryImage);
% Initialize variables to store results
approximatedContours = cell(size(contours));
% Iterate through each contour
for i = 1:length(contours)
    contour = contours{i};
    % Perform polygon approximation using the MPP algorithm
    approximatedContour = mpp(contour);
    % Store results
    approximatedContours{i} = approximatedContour;
end
% Display both images side by side
figure;
% Original binary image
subplot(1, 2, 1);
imshow(binaryImage);
title('Original Binary Image');
hold on;
for i = 0:imageSize(2)
    line([i+0.5, i+0.5], [0.5, imageSize(1)+0.5], 'Color', [1 1 1], 'LineWidth', 2);
end
for i = 0:imageSize(1)
    line([0.5, imageSize(2)+0.5], [i+0.5, i+0.5], 'Color', [1 1 1], 'LineWidth', 2);
end
hold off;
% Display the approximated polygon using the MPP algorithm
subplot(1, 2, 2);
imshow(binaryImage); % Original binary image as the background
hold on;
for i = 0:imageSize(2)
    line([i+0.5, i+0.5], [0.5, imageSize(1)+0.5], 'Color', [1 1 1], 'LineWidth', 2);
end
for i = 0:imageSize(1)
    line([0.5, imageSize(2)+0.5], [i+0.5, i+0.5], 'Color', [1 1 1], 'LineWidth', 2);
end

for i = 1:length(approximatedContours)
    polygon = approximatedContours{i};
    plot(polygon(:, 2), polygon(:, 1), 'r', 'LineWidth', 2);
end

hold off;
title('Minimum Perimeter Polygon(MPP) Approx.');

% MPP algorithm for polygon approximation
function simplified = mpp(points)
    % Initialize simplified polygon
    simplified = [];
    
    % Iterate through the points to determine vertices
    for i = 1:size(points, 1)
        % Check if the current point is a valid vertex
        if isValidVertex(points, i)
            simplified = [simplified; points(i, :)];
        end
    end
end

% Check if a point is a valid vertex based on MPP rules
function valid = isValidVertex(points, index)
    % Check if the point satisfies the MPP rules for valid vertices
    valid = false;
    
    % Get the coordinates of the current, previous, and next points
    currentPoint = points(index, :);
    prevIndex = mod(index - 2, size(points, 1)) + 1;
    nextIndex = mod(index, size(points, 1)) + 1;
    prevPoint = points(prevIndex, :);
    nextPoint = points(nextIndex, :);
    
    % Check if the current point is a white (convex) vertex or a black (mirrored concave) vertex
    isWhiteCurrent = isWhite(currentPoint, prevPoint, nextPoint);
    isBlackCurrent = isBlack(currentPoint, prevPoint, nextPoint);
    
    % Validate the point based on MPP rules
    if isWhiteCurrent || isBlackCurrent
        valid = true;
    end
end

% Check if a point is a white (W) vertex (convex)
function result = isWhite(current, previous, next)
    result = true;
    if ~isConvex(current, previous, next)
        result = false;
    end
end

% Check if a point is a black (B) vertex (mirrored concave)
function result = isBlack(current, previous, next)
    result = true;
    if ~isConcave(current, previous, next)
        result = false;
    end
end

% Check if a point is convex based on the cross product of vectors
function result = isConvex(current, previous, next)
    vector1 = previous - current;
    vector2 = next - current;
    crossProduct = vector1(1) * vector2(2) - vector1(2) * vector2(1);
    angle = atan2d(crossProduct, dot(vector1, vector2));
    result = angle < 180;
end

function result = isConcave(current, previous, next)
    vector1 = previous - current;
    vector2 = next - current;
    crossProduct = vector1(1) * vector2(2) - vector1(2) * vector2(1);
    angle = atan2d(crossProduct, dot(vector1, vector2));
    result = angle >= 180;
end
