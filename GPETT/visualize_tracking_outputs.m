function [  ] = visualize_tracking_outputs( estState, estQuat, estStateVar, measArray...
    , groundTruth, time, basisAngleArray )
% VISUALIZE_TRACKING_OUTPUTS visualizes the 3D extent estimate with the
% associated uncertainty information. Besides, if the ground truth is
% available, it is plotted as well.

% Determine whether the ground truth information is provided
if isempty(groundTruth)
    isGTAvailable = 0;
else
    isGTAvailable = 1;
end

% Extract relevant information from the updated state
estCenter = estState(1:3);
% estVelocity = estState(4:6);
estExtent = estState(13:end);

% For debugging
% estExtent(estExtent<0) = 0;

% Extract extent variances
stdState = sqrt(estStateVar);
stdExtent = stdState(13:end);

% Extract ground truth values of the position and quaternions
if isGTAvailable
    objType = groundTruth.objectDescription(1);
    objParam = groundTruth.objectDescription(2:end);
    gtKinematics = groundTruth.dataLog(abs(groundTruth.dataLog(:,1)-time)<1e-10, 2:end);
    gtCenter =gtKinematics(1:3)';
    gtQuat = gtKinematics(4:7)';
end

cla(); % Clear the current axis

% Plot the point measurements
plot3(measArray(:,1), measArray(:,2), measArray(:,3), 'r+', 'LineWidth', 2, 'MarkerSize', 5);

% Plot the extent estimates with their uncertainty
plot_extent_with_uncertainty(estCenter, estQuat, estExtent, stdExtent, basisAngleArray);

% Plot the ground truth
if isGTAvailable
    switch objType
        case 1 % Sphere
            gTPlotHandle = plot_sphere(gtCenter, objParam, [1 0 1]);
        case 2 % Box
            gTPlotHandle = plot_box(gtCenter, gtQuat, objParam, [1 0 1]);
    end
end

% Attach a legend
legend('Measurements', 'Estimated surface', '1-sigma outer surface', 'Ground truth');

% Arrange lighting in the environment
delete(findall(gcf,'Type','light'))
lightangle(-60,60)
lighting phong;
material shiny;
drawnow();
end

function plotHandles =  plot_extent_with_uncertainty (objCenter, quaternion, extent, stdExtent, basisAngles)

% Calculate the rotation matrix from Local to Global by the quternions
R_from_G_to_L = rotation_matrix_from_global_to_local(quaternion);
R_from_L_to_G = transpose(R_from_G_to_L);

% Convert the extent into cartesion coordinates in Global frame
extentCart_G = transform_extent_from_local_to_global(extent, basisAngles...
    , objCenter, R_from_L_to_G);

% Extract x,y,z coordinates of extent estimates
xExtent = transpose(extentCart_G(1, :));
yExtent = transpose(extentCart_G(2, :));
zExtent = transpose(extentCart_G(3, :));

%% Plot the estimated surface
% First, construct the connectivity list considering the unit sphere
% evaluated at the basis angles
sphere = transform_extent_from_local_to_global(ones(length(basisAngles),1), basisAngles...
    , objCenter, R_from_L_to_G);
connectiviyList = boundary(sphere(1,:)', sphere(2,:)', sphere(3,:)', 0.2); % Perform triangulation

% Then, plot the estimated surface by triangulation
% sEst = trisurf(connectiviyList, xExtent, yExtent, zExtent, 'FaceAlpha', 0.6, 'FaceColor', [0.3 0.75 0.93]); % Plot the triangles
triangles = boundary(extentCart_G', 0.5);
sEst = trisurf(triangles, xExtent, yExtent, zExtent, 'FaceAlpha', 0.6, 'FaceColor', [0.3 0.75 0.93]); % Plot the triangles
sEst.EdgeColor = 'none';

%% Plot outer and inner surfaces regarding the uncertainty of estimates
extentInner = extent - sign(extent) .* min(abs(extent), stdExtent);
extentOuter = extent + sign(extent) .* stdExtent;

% Convert inner extent into cartesion coordinates in Global frame
innerCart_G = transform_extent_from_local_to_global(extentInner, basisAngles...
    , objCenter, R_from_L_to_G);

% Convert outer extent into cartesion coordinates in Global frame
outerCart_G =  transform_extent_from_local_to_global(extentOuter, basisAngles...
    , objCenter, R_from_L_to_G);

% Plot the inner surface
% sInner = trisurf(connectiviyList, innerCart_G(1,:)', innerCart_G(2,:)', innerCart_G(3,:)',...
%     'FaceAlpha', 0.2, 'FaceColor', 'r'); % Plot the triangles
% sInner.EdgeColor = 'none';

% Plot the outer surface
triangles = boundary(outerCart_G', 0.5);
sOuter = trisurf(triangles, outerCart_G(1,:)', outerCart_G(2,:)', outerCart_G(3,:)',...
    'FaceAlpha', 0.1, 'FaceColor', [1 0.9 0]); % Plot the triangles
sOuter.EdgeColor = 'none';

% Return the plot handles
plotHandles = [sEst, sOuter];

    function  extentGlobal = transform_extent_from_local_to_global(extentLocal, basisAngles, objectCenter, rotMatrix)
        
        % Convert the extent into cartesion coordinates in Local frame
        [x_L, y_L, z_L] = sph2cart(basisAngles(:,1), basisAngles(:,2), extentLocal);
        cart_L = [x_L'; y_L'; z_L'];
        
        % Transform the extent to the Global frame
        extentGlobal = objectCenter + rotMatrix * cart_L;
    end
end

function plotHandles = plot_box(origin, quaternions, objectParameters, color)
% CUBE_PLOT plots a cube
% INPUTS:
%   origin = set origin point for the cube in the form of [x,y,z]
%   edgeLength = length of an edge
%   quaternions = [q0 q1 q2 q3] to determine orientation wrt global frame
%   color  = STRING, the color patched for the cube.
%         List of colors
%         b blue
%         g green
%         r red
%         c cyan
%         m magenta
%         y yellow
%         k black
%         w white

% Extract object parameters
length = objectParameters(1);
height = objectParameters(2);
width = objectParameters(3);

% Define the vertices of the box
vertices_L = [length width -height;
    -length width -height;
    -length width height;
    length width height;
    -length -width height;
    length -width height;
    length -width -height;
    -length -width -height]* 0.5;

R_from_G_to_L = rotation_matrix_from_global_to_local(quaternions);
R_from_L_to_G = transpose(R_from_G_to_L);
vertices_G = transpose(R_from_L_to_G * vertices_L');

%  Define the faces of the box
fac = [1 2 3 4;
    4 3 5 6;
    6 7 8 5;
    1 2 8 7;
    6 7 1 4;
    2 3 5 8];
box = [vertices_G(:,1)+origin(1),vertices_G(:,2)+origin(2),vertices_G(:,3)+origin(3)];
% plotHandles = patch('Faces',fac,'Vertices',box,'FaceColor',color, 'EdgeColor'...
%     , 'none', 'LineWidth', 1,'LineStyle', ':','FaceAlpha', 0.35);

plotHandles = patch('Faces',fac,'Vertices',box,'FaceColor',color, 'EdgeColor'...
    , [0.3 0.3 0.3], 'LineWidth', 1,'FaceAlpha', 0);
end

function plotHandles = plot_sphere(center, objectParameters, color)
% Plots a sphere

% Extract object parameters
radius = objectParameters;% Ellipsoid

% Produce cartesian points on the unit sphere
[xSphere, ySphere, zSphere] = sphere(128);

% Scale and translate these points
xSphere = xSphere*radius + center(1);
ySphere = ySphere*radius + center(2);
zSphere = zSphere*radius + center(3);

% Plot the surface
h = surfl(xSphere, ySphere, zSphere);
h.FaceColor = color;
h.FaceAlpha = 0.3;
h.EdgeColor = 'none';

% Return the surface handle
plotHandles = h;
end



