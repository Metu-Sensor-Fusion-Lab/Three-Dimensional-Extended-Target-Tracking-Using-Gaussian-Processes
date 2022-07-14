function [recontructedObject3D] = reconstruct_from_projections(center3D, quat... 
    , extent1, extent2, extent3, basisAngleArray)
% This function reconstructs the 3D object from its three projected contours. 
% The 3D object is represented by 3D points. 

% Contour1 is on XY-plane, Contour2 is on XZ-plane, Contour3 is on YZ-plane of the local frame
P1 = [1 0 0; 0 1 0];    % Projects to Plane-1 (X-Y plane)
P2 = [1 0 0; 0 0 1];    % Projects to Plane-2 (X-Z plane)
P3 = [0 1 0; 0 0 1];    % Projects to Plane-3 (Y-Z plane)
numBasisAngles = size(basisAngleArray, 1);
margin = 0.5;           % Extention from the contour limits (in m)
numLayers = 20;         % Number of test layers (the contours will be repeated for this many times)

% Determining Contour1 in Cartesian coordinates
[x1, y1] = pol2cart(basisAngleArray, extent1);

% Determining Contour2 in Cartesian coordinates
[x2, z2] = pol2cart(basisAngleArray, extent2);

% Determining Contour3 in Cartesian coordinates
[y3, z3] = pol2cart(basisAngleArray, extent3);

% Construct the contours in Cartesian coordinates
contour1 = [x1 y1; [x1(1) y1(1)]];  % Repeat the first vertex to close the loop
contour2 = [x2 z2; [x2(1) z2(1)]];
contour3 = [y3 z3; [y3(1) z3(1)]];

% Determine maximums and minimums from the 2D contours
xmin = min([x1; x2]) - margin;
xmax = max([x1; x2]) + margin;
ymin = min([y1; y3]) - margin;
ymax = max([y1; y3]) + margin;
zmin = min([z2; z3]) - margin;
zmax = max([z2; z3]) + margin;

% Build a linspaced array for each axis
xArray = linspace(xmin, xmax, numLayers);
xArray = repmat(xArray, numBasisAngles, 1);
xArray = xArray(:);
yArray = linspace(ymin, ymax, numLayers);
yArray = repmat(yArray, numBasisAngles, 1);
yArray = yArray(:);
zArray = linspace(zmin, zmax, numLayers);
zArray = repmat(zArray, numBasisAngles, 1);
zArray = zArray(:);

% Build a rough shape by repeated contours
contour1Repeated = [repmat(x1, numLayers, 1) repmat(y1, numLayers, 1) zArray];
contour2Repeated = [repmat(x2, numLayers, 1) yArray repmat(z2, numLayers, 1)];
contour3Repeated = [xArray repmat(y3, numLayers, 1) repmat(z3, numLayers, 1)];

testPoints = [contour1Repeated; contour2Repeated; contour3Repeated];

%% Test points are to be eleminated by the following carving process
% Plane 1
testPointsProjected = P1 * testPoints';   % Project all points on Plane1
testPointsProjected = testPointsProjected';
% Delete the points that are outside Contour1
in = inpolygon(testPointsProjected(:, 1), testPointsProjected(:, 2), contour1(:,1), contour1(:,2)); 
testPoints(~in, :) = []; 

% Plane 2
testPointsProjected = P2 * testPoints';   % Project all points on Plane2
testPointsProjected = testPointsProjected';
% Delete the points that are outside Contour2
in = inpolygon(testPointsProjected(:, 1), testPointsProjected(:, 2), contour2(:,1), contour2(:,2)); 
testPoints(~in, :) = []; 

% Plane3
testPointsProjected = P3 * testPoints';
testPointsProjected = testPointsProjected';
% Delete the points that are outside Contour3
in = inpolygon(testPointsProjected(:, 1), testPointsProjected(:, 2), contour3(:,1), contour3(:,2)); 
testPoints(~in, :) = []; 

% Rotate the points regarding the orientation of the 3D object
R_from_G_to_L = rotation_matrix_from_global_to_local(quat);
R_from_L_to_G = transpose(R_from_G_to_L);
testPoints_G = transpose(R_from_L_to_G * testPoints');

% Shift the points to the object center
recontructedObject3D = testPoints_G + center3D'; 

end