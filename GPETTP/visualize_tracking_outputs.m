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
numBasisAngles = size(basisAngleArray,1);
estExtent_P1 = estState(13:12+numBasisAngles);
estExtent_P2 = estState(13+numBasisAngles:12+2*numBasisAngles);
estExtent_P3 = estState(13+2*numBasisAngles:12+3*numBasisAngles);

% Extract extent variances
stdState = sqrt(estStateVar);
stdExtent_P1 = stdState(13:12+numBasisAngles);
stdExtent_P2 = stdState(13+numBasisAngles:12+2*numBasisAngles);
stdExtent_P3 = stdState(13+2*numBasisAngles:12+3*numBasisAngles);

% Extract ground truth values of the position and quaternions
if isGTAvailable
    objType = groundTruth.objectDescription(1);
    objParam = groundTruth.objectDescription(2:end);
    gtKinematics = groundTruth.dataLog(abs(groundTruth.dataLog(:,1)-time)<1e-10, 2:end);
    gtCenter =gtKinematics(1:3)';
    gtQuat = gtKinematics(4:7)';
end

% Contruct the 3D object from the projection contours
reconstructedObject3D = reconstruct_from_projections(estCenter, estQuat...
    , estExtent_P1, estExtent_P2, estExtent_P3, basisAngleArray);

cla(); % clear the current axis

% Plot the point measurements
plot3(measArray(:,1), measArray(:,2), measArray(:,3), 'r+', 'LineWidth', 2, 'MarkerSize', 5);

% Plot reconstructed object
plot_reconstructed_3D_object(reconstructedObject3D);

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
legend('Measurements', 'Estimated surface', 'Ground truth');

% Arrange lighting in the environment
delete(findall(gcf,'Type','light'))
lightangle(-60,60)
lighting phong;
material shiny;
drawnow();
end

function plotHandles = plot_box(origin, quaternions, objectParameters, color)

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

plotHandles = patch('Faces',fac,'Vertices',box,'FaceColor',color, 'EdgeColor'...
    , [0.3 0.3 0.3], 'LineWidth', 1,'FaceAlpha', 0);
end

function plotHandles = plot_sphere(center, objectParameters, color)
% Plots a sphere

% Extract object parameters
radius = objectParameters;

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

function [ plotHandles ] = plot_reconstructed_3D_object( object3DPoints )

triangles = boundary(object3DPoints, 0.5);                                      % perform triangulation
s = trisurf(triangles, object3DPoints(:,1), object3DPoints(:,2)...
    , object3DPoints(:,3), 'FaceAlpha', 0.6, 'FaceColor', [0.3 0.75 0.93]);     % plot the triangulated surface
s.EdgeColor = 'none';

plotHandles = s;

end
