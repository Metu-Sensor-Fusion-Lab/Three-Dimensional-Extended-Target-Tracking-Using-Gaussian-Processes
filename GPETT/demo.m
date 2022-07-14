% This file presents a demonstration of the extended tracker proposed in 
% Section V of the following study.

% For the details, please see:
%       M. Kumru and E. Özkan, “Three-Dimensional Extended Object Tracking 
%       and Shape Learning Using Gaussian Processes”, 
%       IEEE Transactions on Aerospace and Electronic Systems, 2021.
%       https://ieeexplore.ieee.org/abstract/document/9382878/

% Dependencies:
%       "Spheretri" package by Peter Gagarinov
%       It is utilized to produce uniformly distributed points on the unit sphere.
%       https://www.mathworks.com/matlabcentral/fileexchange/58453-spheretri

clc;
clear;
close all;


%% Paramaters
% Simulation Parameters
expDuration = 30;   % (in seconds)
T = 0.1;            % sampling time (in seconds)
eps = 1e-6;         % a small scalar

% Measurement Parameters
numMeasPerInstant = 20;
paramMeas = {numMeasPerInstant, expDuration, T};

% Gaussian Process (GP) Parameters
minNumBasisAngles = 200;        % minimum number of basis angles on which the extent is maintained
meanGP = 0;                     % mean of the GP
stdPriorGP = 1;                 % prior standard deviation of the GP
stdRadiusGP = 0.2;              % standard deviation of the radius
scaleLengthGP = pi/8;           % lengthscale
kernelTypeGP = 1;               % 1:Squared exponential , 2: Matern3_2, 3: Matern5_2
isObjSymmetric = 0;             % it is used to adjust the kernel function accordingly
stdMeasGP = 0.1;                % standard deviation of the measurements (used in the GP model)
paramGP = {kernelTypeGP, stdPriorGP, stdRadiusGP, scaleLengthGP, isObjSymmetric, stdMeasGP, meanGP};

% EKF Paramaters
stdCenter = 1e-1;               % std dev of the process noise for object center
stdAngVel = 1e-1;               % std dev of the process noise for angular velocity
lambda = 0.99;                   

% Determine the Basis Angles
[basisVertices, ~] = spheretri(minNumBasisAngles);      % produces evenly spaced points on the sphere
[azimuthBasisArray, elevationBasisArray, ~] = cart2sph(basisVertices(:,1), basisVertices(:,2)...
    , basisVertices(:,3));
azimuthBasisArray = mod(azimuthBasisArray, 2*pi);       % to make it consistent with the convention
numBasisAngles = size(basisVertices, 1);
% Arrange the basis array so that azimuth and elevation are in ascending order
basisAngleArray = [azimuthBasisArray elevationBasisArray];
basisAngleArray = sortrows(basisAngleArray, 2);
basisAngleArray = sortrows(basisAngleArray, 1);


%% Load Measurements and Ground Truth 
[measComplete, groundTruth] = get_measurements(paramMeas);


%% Specify the Initial Distribution
% Initial covariance
P0_c = 0.1 * eye(3);
P0_v = 0.1 * eye(3);
P0_a = 1e-5 * eye(3,3);
P0_w = 1e-0 * eye(3,3);
P0_extent = compute_GP_covariance(basisAngleArray, basisAngleArray, paramGP);  
P0 = blkdiag(P0_c, P0_v, P0_a, P0_w, P0_extent);

% Initial mean
c0 = groundTruth.dataLog(1, 2:4)'+ sqrt(0.1)*randn(3,1);        % initial linear velocity
v0 = [0.5 0 0]'+ sqrt(0.1)*randn(3,1);                          % initial linear velocity
q0 = [0 0 0 1]';                                                % initial quaternions 
a0 = [0 0 0]';                                                  % initial orientation deviation
w0 = [0 0 0]';                                                  % initial angular velocity
f0 = meanGP * ones(numBasisAngles, 1);                          % initial extent (is initialized regarding the GP model)

% Initialize the filter estimate
estState = [c0; v0; a0; w0; f0];
estQuat = q0;
estStateCov = P0;


%% Define the Process Model
F_lin = kron([1 T; 0 1], eye(3));       % constant velocity model for linear motion
F_rot = eye(6,6);                       % a substitute matrix for the rotational dynamics (it will be computed in the filter at each recursion)
F_f = eye(numBasisAngles);      
F = blkdiag(F_lin, F_rot, F_f);

Q_lin = kron([T^3/3 T^2/2; T^2/2 T], stdCenter^2*eye(3));   % process covariance of linear motion
Q_rot = eye(6,6);                                           % a substitute matrix for the covariance of the rotational motion (it will be computed in the filter at each recursion)
Q_extent = zeros(numBasisAngles);                           % predicted covariance of the extent will be computed within the filter according to Eqn. 20 
Q = blkdiag(Q_lin, Q_rot, Q_extent);

processModel = {F, Q, lambda};

 
%% Open a Figure (to illustrate the tracking outputs)
figure;
grid on; hold on; axis equal;
view(45, 25);
title('Extended Target Tracking Results');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');


%% GPETT3D LOOP
time = 0;
numInstants = ceil(expDuration/ T);
filterParam = {processModel, paramGP, basisAngleArray, T, stdAngVel, eps}; 

for k = 1:numInstants
    
    % Extract current measurements
    curMeasArray = measComplete(abs(measComplete(:, 1)-time)<eps, 2:4);    
    
    % Call the GPETT3D filter
    [estState, estQuat, estStateCov] = filter_GPETT3D(estState, estQuat...
        ,estStateCov, curMeasArray, filterParam);
    
    % Visualize the results
    visualize_tracking_outputs(estState, estQuat, diag(estStateCov), curMeasArray, groundTruth...
        ,time, basisAngleArray);        

    % Update time
    time = time + T;
end