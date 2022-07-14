% This file presents a demonstration of the extended tracker proposed in 
% Section VII of the following study.

% For the details, please see:
%       M. Kumru and E. Özkan, “Three-Dimensional Extended Object Tracking 
%       and Shape Learning Using Gaussian Processes”, 
%       IEEE Transactions on Aerospace and Electronic Systems, 2021.
%       https://ieeexplore.ieee.org/abstract/document/9382878/

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
numBasisAngles = 50;            % number of basis angles on which extent of each projection is maintained
meanGP = 0;                     % mean of the GP
stdPriorGP = 1;                 % prior standard deviation of the GP
stdRadiusGP = 0.2;              % standard deviation of the radius
scaleLengthGP = pi/5;           % scale length
kernelTypeGP = 1;               % 1:Squared exponential, 2: Matern3_2, 3: Matern5_2
isObjSymmetric = 0;             % it is used to adjust the kernel function accordingly
stdMeasGP = 0.2;                % standard deviation of the measurements (used in the GP model)
paramGP1 = {kernelTypeGP, stdPriorGP, stdRadiusGP, scaleLengthGP, 0, stdMeasGP, meanGP};
paramGP2 = {kernelTypeGP, stdPriorGP, stdRadiusGP, scaleLengthGP, 0, stdMeasGP, meanGP};
paramGP3 = {kernelTypeGP, stdPriorGP, stdRadiusGP, scaleLengthGP, 0, stdMeasGP, meanGP};
paramGP = {paramGP1, paramGP2, paramGP3};

% EKF Paramaters
stdCenter = 1e-1;               % std dev of the object center
stdAngVel = 4e-1;               % std dev of each Euler angle
lambda = 0.99;

% Parameters of the Scaling
meanS = 5/6;                    % mean of the scaling parameter
stdS = sqrt(1/18);              % standard deviation of the scaling parameter
paramScaling = [meanS stdS];

% Determine the Basis Angles for 2D Contour
basisAngleArray = transpose(linspace(0, 2*pi, numBasisAngles+1));
basisAngleArray(end) = [];      % Trash the end point to eleminate repetition at 2*pi


%% Load Measurements and Ground Truth
[measComplete, groundTruth] = get_measurements(paramMeas);


%% Determine the Initial Distribution
% Initial covariance
P0_c = 0.1 * eye(3);
P0_v = 0.1 * eye(3);
P0_a = 1e-4 * eye(3,3);
P0_w = 1e-4 * eye(3,3);
P0_extent1 = compute_GP_covariance(basisAngleArray, basisAngleArray, paramGP{1});
P0_extent2 = compute_GP_covariance(basisAngleArray, basisAngleArray, paramGP{2});
P0_extent3 = compute_GP_covariance(basisAngleArray, basisAngleArray, paramGP{3});
P0 = blkdiag(P0_c, P0_v,  P0_a, P0_w, P0_extent1, P0_extent2, P0_extent3);

% Initial mean
c0 = groundTruth.dataLog(1, 2:4)' + sqrt(0.1)*randn(3,1);       % initial position
v0 = [0.5 0 0]' + sqrt(0.1)*randn(3,1);                         % initial linear velocity
q0 = [0 0 0 1]';                                                % initial quaternion
a0 = [0 0 0]';                                                  % initial orientation deviation vector
w0 = [0 0 0]';                                                  % initial angular  velocity
f1_0 = meanGP * ones(numBasisAngles, 1);                        % initial extent for the first projection (is initialized regarding the GP model)
f2_0 = f1_0;                                                    % initial extent for the second projection
f3_0 = f1_0;                                                    % initial extent for the third projection
f0 = [f1_0; f2_0; f3_0];                                        % initial value of the state vector

% Initialize the filter estimate
estState = [c0; v0; a0; w0; f0];
estQuat = q0;
estStateCov = P0;


%% Define the Process Model
F_lin = kron([1 T; 0 1], eye(3,3));         % constant velocity model for linear motion
F_rot = eye(6,6);                           % a substitute matrix for the rotational dynamics (it will be computed in the filter at each recursion)
F_f1 = eye(numBasisAngles);      
F_f2 = F_f1;
F_f3 = F_f1;
F_f = blkdiag(F_f1, F_f2, F_f3);
F = blkdiag(F_lin, F_rot, F_f);

Q_lin = kron([T^3/3 T^2/2; T^2/2 T], stdCenter^2*eye(3));    % process noise covariance of linear motion
Q_rot = eye(6,6);                                            % a substitute matrix for the covariance of the rotational motion (it will be computed in the filter at each recursion)
Q_extent1 = zeros(numBasisAngles);                       
Q_extent2 = Q_extent1;
Q_extent3 = Q_extent1;
Q_extent = blkdiag(Q_extent1, Q_extent2, Q_extent3);
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


%% GPETT3DProjection LOOP
time = 0;
numInstants = ceil(expDuration/ T);
filterParam = {processModel, paramGP, paramScaling, basisAngleArray, T, stdAngVel, eps}; 

for k = 1:numInstants
    
    % Extract current measurements
    curMeasArray = measComplete(abs(measComplete(:, 1)-time)<eps, 2:4);
    
    % Call the GPETT3DProjection filter
    [estState, estQuat, estStateCov] = filter_GPETT3DProjection(estState, estQuat...
        , estStateCov, curMeasArray, filterParam);
    
    % Visualize the results
    visualize_tracking_outputs(estState, estQuat, diag(estStateCov), curMeasArray, groundTruth...
        , time, basisAngleArray);
    
    % Update time
    time = time + T;
end
