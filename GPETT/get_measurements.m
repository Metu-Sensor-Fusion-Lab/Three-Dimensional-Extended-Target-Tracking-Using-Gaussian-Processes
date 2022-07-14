function [ measComplete, groundTruth ] = get_measurements( paramMeas )
%GET_MEASUREMENTS provides point measurements. 
% Outputs:
%	MeasComplete: is the complete set of point measurements. Each row is in
%       the form of [timeStamp x y z]
%   groundTruth: is the ground truth information of the corresponding scenario.
%       It is a structure with the following fields: 'objectDescription', 'dataRecord'.
%       'objectDescription' includes type (e.g., sphere, box...) and the parameters
%       of the geometry.
%       Each row of the 'dataRecord' is in the form of [timeStamp center' quaternions']
%       where centerPosition = [cX cY cZ]' and quaternions = [q0 q1 q2 q3]'.

% Input:
%   paramMeas: measurement parameters. It is cell variable including the
%       following information {measSource, numMeasPerInstant, expDuration, T}

% Extract the parameters
numMeasPerInstant = paramMeas{1};
expDuration = paramMeas{2};
T = paramMeas{3};
numInstants = ceil(expDuration/ T);
stdMeas = 0.1;              % Measurement noise std

% Object parameters
objType = 2; % objType = 1:Sphere, 2:Box
switch objType
    case 1 % Sphere
        radius = 2;
        objParameters = radius;
    case 2 % Box
        l = 3;      % Length
        h = 3;      % Height
        w = 3;      % Width
        objParameters = [l h w];
end

% Generate the measurements in Matlab
[measComplete, gtDataLog] = generate_measurements_Matlab...
    (stdMeas, numMeasPerInstant, objType, objParameters, numInstants, T);

groundTruth = struct('objectDescription', [objType, objParameters]...
    , 'dataLog', gtDataLog);

end

function [surfaceMeasurements, groundTruth] = generate_measurements_Matlab...
    (stdMeas, numMeasPerInstant, objType, objParameters, numInstants, T)
% This function produces surface measurements from a 3D object.
% Outputs:
%       surfaceMeasurements = [timeStamp xPosition yPosition zPosition]
%       groundTruth = [timeStamp centerPosition' quaternions']
%           where centerPosition = [cX cY cZ] and quaternions = [q0 q1 q2 q3]
% Inputs:
%       stdMeas: standart deviation of a 3D point measurement on each axis
%       numInstants: number of instants included in the simulation
%       objType: enumerated variable, 1:Sphere, 2:Cube, 3:Ellipsoid,...
%       objParameters: depends on the object, i.e., for Cube it denotes edgeLength
%       experimentDuration: in seconds
%       T: time step between successive siulation instants (in seconds)

% Comment out one of the test scenarios below (A or B)
% A. Straight-line pattern:
% initialQuaternion = [0 ; 0; 0; 1];
% initialPosition = [0; 0; 0];
% linearVelocityMagnitude = 10; % in m/sec
% angularVelocityArray = zeros(numInstants, 3);
% angularVelocityPhi = 0;         % in rad/sec, around X axis of the local frame
% angularVelocityTheta = 0;   % in rad/sec, around Y axis of the local frame
% angularVelocityPsi = 0;    % in rad/sec, around Z axis of the local frame
% angularVelocityArray(:, 1) = angularVelocityPhi; 
% angularVelocityArray(:, 2) = angularVelocityTheta; 
% angularVelocityArray(31:230, 3) = angularVelocityPsi; 

% B. U-turn pattern:
initialQuaternion = [0; 0; 0; 1];
initialPosition = [0; 0; 0];
linearVelocityMagnitude = 0.5;          % in m/sec
angularVelocityArray = zeros(numInstants, 3);
angularVelocityPhi = pi/10;             % in rad/sec, around X axis of the local frame
angularVelocityTheta = pi/15;           % in rad/sec, around Y axis of the local frame
angularVelocityPsi = pi/20;             % in rad/sec, around Z axis of the local frame
angularVelocityArray(30:100, 1) = -angularVelocityPhi; 
angularVelocityArray(170:220, 1) = 0.7*angularVelocityPhi; 
angularVelocityArray(50:120, 2) = angularVelocityTheta; 
angularVelocityArray(190:250, 2) = -1*angularVelocityTheta; 
angularVelocityArray(80:150, 3) = angularVelocityPsi; 
angularVelocityArray(220:260, 3) = 0*angularVelocityPsi; 

%% Simulate the Motion
% Note that only kinematic aspects of motion is considered. No dynamic behaviour is introduced.
posArray = zeros(numInstants, 3); % : [xArray yArray zArray]
quatArray = zeros(numInstants, 4); % : [q0Array q1Array q2Array q3Array]
velArray = zeros(numInstants, 3); 
groundTruth = zeros(numInstants, 14); % : [timeStamp centerPosition' quaternions' angularRate' linearVel']

pos_k = initialPosition; % Initialize position
quat_k = normc(initialQuaternion); % Initialize quaternions
vel_k = linearVelocityMagnitude * [1; 0; 0];
angRate_k = angularVelocityArray(1,:)';

% Record the initial state
posArray(1, :) = pos_k';
velArray(1, :) = vel_k';
quatArray(1, :) = quat_k';
groundTruth(1, :) = [0 pos_k' quat_k' angRate_k' vel_k'];

for k = 1:(numInstants-1)
    curTime = k*T;
    
    % Obtain transformation matrices
    % First, find the linear velocity rotation matrix from local frame to the global
    R_from_G_to_L = rotation_matrix_from_global_to_local(quat_k);
    R_from_L_to_G = transpose(R_from_G_to_L);
    % Then, find the matrix to transform local angular velocities into quaternion diffential
    T_from_L_to_G = angular_velocity_transformation_matrix_from_local_to_global(quat_k);
    
    % Obtain the linear velocity vector expressed in the local frame
    v_L_k = linearVelocityMagnitude * [1; 0; 0]; % It is assumed to be in X direction
    
    % Extract the current angular velocity vector
    w_L_k = transpose(angularVelocityArray(k, :));
    
    % Convert the velocities into global frame
    posDot_G_k = R_from_L_to_G * v_L_k;
    %     quatDot_G_k = T_from_L_to_G * w_L_k;
    quatDot_G_k = 1/2 * [quat_k(4)*w_L_k - cross(w_L_k, quat_k(1:3)); -w_L_k'*quat_k(1:3)];
    
    % Update the position and quaternions
    pos_k = pos_k + T * posDot_G_k;
    quat_k = quat_k + T * quatDot_G_k;
    quat_k = normc(quat_k);
    
    % Record the variables
    posArray(k+1, :) = pos_k';
    velArray(k+1, :) = posDot_G_k';
    quatArray(k+1, :) = quat_k';
    groundTruth(k+1, :) = [curTime pos_k' quat_k' angularVelocityArray(k,:) posDot_G_k'];
end

%% Produce surface measurements
surfaceMeasurements = zeros(numInstants*numMeasPerInstant, 4); % Initialize the measurements
for i = 1:numInstants
    iTime = groundTruth(i, 1);
    iPos = groundTruth(i, 2:4)';
    iQuat = groundTruth(i, 5:8)';
    
    curMeasArray = produce_surface_measurements(iPos, iQuat, numMeasPerInstant, stdMeas, objType,...
        objParameters);
    
    % Insert the measurements to the output matrix with relevant time stamp
    surfaceMeasurements((i-1)*numMeasPerInstant + 1 : (i*numMeasPerInstant), :) = ...
        [ones(numMeasPerInstant, 1)*iTime curMeasArray];
end
end

function [surfaceMeas] = produce_surface_measurements(pos, quat, numMeasurements, stdMeas, objType...
    , objParameters)
% PRODUCE_SURFACE_MEASUREMENTS generates 3D point measurements from the
% surface of an object with specified geometry.
switch objType
    case 1 % Sphere
        radius = objParameters;
        
        % Assume x,y and z are independent Gaussian random variables
        sampleArray = randn(3, numMeasurements);
        % Project each sample onto the sphere by normalizing the length
        sampleArray = normc(sampleArray);
        % Multiply the normalized sample to obtain desired radius
        sampleArray = radius .* sampleArray;
        
        xTrue_L = sampleArray(1,:)';
        yTrue_L = sampleArray(2,:)';
        zTrue_L = sampleArray(3,:)';
        
    case 2 % Box
        l = objParameters(1);
        h = objParameters(2);
        w = objParameters(3);
        
        % First, produce a random array for the selection of the faces
        faceArray = randi([0 5], numMeasurements, 1);
        
        % Initialize the outputs
        xTrue_L = zeros(numMeasurements, 1);
        yTrue_L = zeros(numMeasurements, 1);
        zTrue_L = zeros(numMeasurements, 1);
        
        % Select random points exactly on the specified box
        pointer = 1;
        for iFace = 0:5
            num = sum(faceArray == iFace);
            switch iFace
                case 0 % on X-Y plane, Z is positive
                    xTrue_L(pointer:(pointer+num-1)) = -l/2 + l*rand(num, 1);
                    yTrue_L(pointer:(pointer+num-1)) = -w/2 + w*rand(num, 1);
                    zTrue_L(pointer:(pointer+num-1)) = h/2*ones(num,1);
                case 1 % on X-Y plane, Z is negative
                    xTrue_L(pointer:(pointer+num-1)) = -l/2 + l*rand(num, 1);
                    yTrue_L(pointer:(pointer+num-1)) = -w/2 + w*rand(num, 1);
                    zTrue_L(pointer:(pointer+num-1)) = -h/2*ones(num,1);
                case 2 % on X-Z plane, Y is positive
                    xTrue_L(pointer:(pointer+num-1)) = -l/2 + l*rand(num, 1);
                    yTrue_L(pointer:(pointer+num-1)) = w/2*ones(num,1);
                    zTrue_L(pointer:(pointer+num-1)) = -h/2 + h*rand(num, 1);
                case 3 % on X-Z plane, Y is negative
                    xTrue_L(pointer:(pointer+num-1)) = -l/2 + l*rand(num, 1);
                    yTrue_L(pointer:(pointer+num-1)) = -w/2*ones(num,1);
                    zTrue_L(pointer:(pointer+num-1)) = -h/2 + h*rand(num, 1);
                case 4 % on Y-Z plane, X is positive
                    xTrue_L(pointer:(pointer+num-1)) = l/2*ones(num,1);
                    yTrue_L(pointer:(pointer+num-1)) = -w/2 + w*rand(num, 1);
                    zTrue_L(pointer:(pointer+num-1)) = -h/2 + h*rand(num, 1);
                case 5 % on Y-Z plane, X is negative
                    xTrue_L(pointer:(pointer+num-1)) = -l/2*ones(num,1);
                    yTrue_L(pointer:(pointer+num-1)) = -w/2 + w*rand(num, 1);
                    zTrue_L(pointer:(pointer+num-1)) = -h/2 + h*rand(num, 1);
            end
            pointer = pointer + num;
        end   
end

% Add Gaussian noise on selected points
xNoisy_L = xTrue_L + stdMeas * randn(numMeasurements, 1);
yNoisy_L = yTrue_L + stdMeas * randn(numMeasurements, 1);
zNoisy_L = zTrue_L + stdMeas * randn(numMeasurements, 1);
meas_L = [xNoisy_L'; yNoisy_L'; zNoisy_L'];

% Obtain the rotation matrix from local frame to the global
R_from_G_to_L = rotation_matrix_from_global_to_local(quat);
R_from_L_to_G = transpose(R_from_G_to_L);

% Rotate 3D points to the global frame
meas_G = transpose(R_from_L_to_G * meas_L);

% Translate the measurements to the center of the object
surfaceMeas = [meas_G(:,1)+pos(1) meas_G(:,2)+pos(2) meas_G(:,3)+pos(3)];
end

function [ transformationMatrix ] = angular_velocity_transformation_matrix_from_local_to_global( quat )
% The function computes the transformation matrix from Local (body) frame to the Global (inertial)
% frame.
% Note that this matrix is to multiply the angular velocity vector expressed in Local from from LEFT,
% i.e., quatDot  =  transformationMatrix  *  wLocal

q0 = quat(1);
q1 = quat(2);
q2 = quat(3);
q3 = quat(4);

% Rotation matrix from global to local R_from_G_to_L
transformationMatrix = 1/2 * [-q1 -q2 -q3; q0 -q3 q2; q3 q0 -q1; -q2 q1 q0];
end
