function [ estState, estQuat, estStateCov ] = filter_GPETT3DProjection( prevEstState...
    , prevEstQuat, prevEstStateCov, measArray, filterParam )

% Extract the parameters 
processModel = filterParam{1};
paramGP = filterParam{2};
paramScaling = filterParam{3};
basisAngleArray = filterParam{4};
T = filterParam{5};
stdAngVel = filterParam{6};
eps = filterParam{7};

F = processModel{1};
Q = processModel{2};
lambda = processModel{3};

% Extract scaling parameters
meanS = paramScaling(1);
stdS = paramScaling(2);

stdMeasGP = paramGP{1}{6};

% Compute the inverse of P0_extent since it will repeatedly be used in the system dynamic model
persistent inv_P0_extent1 inv_P0_extent2 inv_P0_extent3;
if isempty(inv_P0_extent1)
    inv_P0_extent1 = compute_inverse_GP_covariance(basisAngleArray, paramGP{1}, eps);
    inv_P0_extent2 = compute_inverse_GP_covariance(basisAngleArray, paramGP{2}, eps);
    inv_P0_extent3 = compute_inverse_GP_covariance(basisAngleArray, paramGP{3}, eps);
end

% Defining the Projection Matrices
% i.e., let X be a vector in R^3, x_p1 = P1*X is the projection of X on Plane-1
persistent P1 P2 P3;
if isempty(P1)
    P1 = [1 0 0; 0 1 0];    % projects to Plane-1 (X-Y plane)
    P2 = [1 0 0; 0 0 1];    % projects to Plane-2 (X-Z plane)
    P3 = [0 1 0; 0 0 1];    % projects to Plane-3 (Y-Z plane)
end


%% Process Update
% Compute the rotational dynamic model
angVelEst = prevEstState(10:12);
[FRot, QRot] = compute_rotational_dynamics(angVelEst, stdAngVel, T);

% Substitute the rotational model
F(7:12, 7:12) = FRot;
Q(7:12, 7:12) = QRot;

% Compute predicted state and covariance
predState = F * prevEstState;
predQuatDev = predState(7:9);
predStateCov = F * prevEstStateCov * F' + Q;

% Dynamic model for maximum entropy in extent
predStateCov(13:end, 13:end) = prevEstStateCov(13:end, 13:end)/ lambda;

% Extract relevant information from the predicted state
predCenter = predState(1:3);
predVelocity = predState(4:6);
numBasisAngles = size(basisAngleArray,1);
predExtent_P1 = predState(13:12+numBasisAngles);
predExtent_P2 = predState(13+numBasisAngles:12+2*numBasisAngles);
predExtent_P3 = predState(13+2*numBasisAngles:12+3*numBasisAngles);


%% Measurement Update
curNumMeas = size(measArray, 1);
numStates = size(F,1);

% In the below loops, three operations are performed:
% 1. Compute measurement predictions of the predicted state by the original
% nonlinear measurement model.
% 2. Obtain the corresponding measurement covariance matrix.
% 3. Linearize the measurement model to employ in EKF equations.
predMeas_P1 = zeros(curNumMeas * 2, 1);             % of the form [x1; y1;...; xn; yn]
predMeas_P2 = zeros(curNumMeas * 2, 1);
predMeas_P3 = zeros(curNumMeas * 2, 1);
measCov_P1 = zeros(curNumMeas * 2);                 % covariance matrix of projection 1
measCov_P2 = zeros(curNumMeas * 2);
measCov_P3 = zeros(curNumMeas * 2);
dH_dX_P1 = zeros(curNumMeas * 2, numStates);        % partial derivatives wrt states
dH_dX_P2 = zeros(curNumMeas * 2, numStates);
dH_dX_P3 = zeros(curNumMeas * 2, numStates);
dH_dY_P1 = zeros(curNumMeas * 2, 3);                % partial derivatives wrt 3D measurement
dH_dY_P2 = zeros(curNumMeas * 2, 3);                % partial derivatives wrt 3D measurement
dH_dY_P3 = zeros(curNumMeas * 2, 3);                % partial derivatives wrt 3D measurement

% Compute the rotation matrix. It will be used to transform the
% measurements into the local frame
R_from_G_to_L = rotation_matrix_from_global_to_local(apply_quat_deviation(predQuatDev, prevEstQuat));

% Plane 1
for i = 1:curNumMeas
    iMeas3D_G = transpose(measArray(i, :));
    
    % Compute the following variables to exploit in measurement model
    [p, H_f, covMeasBasis, angleMeas_L] = compute_measurement_model_components...
        (iMeas3D_G, predState, prevEstQuat, P1, paramGP{1}...
        , basisAngleArray, inv_P0_extent1);
    H_tilda = p * H_f;
    
    % Obtain the prediction for the nonlinear implicit measurement
    % model, h(x,y) = 0
    iMeasModelPredicted = P1 * R_from_G_to_L * (iMeas3D_G - predCenter)...
        - meanS * H_tilda * predExtent_P1;
    predMeas_P1(1+(i-1)*2 : i*2) = iMeasModelPredicted;
    
    % Obtain the covariance of the current measurement
    % Definition: iCovMeas = k(uk,uk), uk: argument of the current measurement
    iCovMeas = compute_GP_covariance(angleMeas_L, angleMeas_L, paramGP{1});
    % Definition: iR_f = k(uk,uk) - k(uk, uf) * inv(K(uf,uf)) * k(uf, uk)
    iR_f = iCovMeas - covMeasBasis * inv_P0_extent1 * covMeasBasis';
    % Compute the total covariance of this specific measurement
    iR =  stdMeasGP^2 * eye(2) + stdS^2*p*iR_f*p' + stdS^2*H_tilda...
        *(predExtent_P1*predExtent_P1')* H_tilda';
    % Plug this matrix into the overall covariance matrix
    measCov_P1(((i-1)*2+1):i*2, ((i-1)*2+1):i*2) = iR;
    
    %% Obtain the partial derivatives of the measurement model
    [H_center, H_velocity, H_quatDev, H_angVel, H_extent, H_meas] = compute_linearized_measurement_matrix...
        (iMeas3D_G, predState, prevEstQuat, P1, predExtent_P1...
        , paramGP{1}, basisAngleArray, inv_P0_extent1, meanS, eps);
    
    dH_dX_P1(1+(i-1)*2:i*2, :) = [H_center H_velocity H_quatDev H_angVel...
        H_extent zeros(2, numBasisAngles)  zeros(2, numBasisAngles)];
    
    dH_dY_P1(1+(i-1)*2:i*2, :) = H_meas;
end

% Plane 2
for i = 1:curNumMeas
    iMeas3D_G = transpose(measArray(i, :));
    
    % Compute the following variables to exploit in measurement model
    [p, H_f, covMeasBasis, angleMeas_L] = compute_measurement_model_components...
        (iMeas3D_G, predState, prevEstQuat, P2, paramGP{2}...
        , basisAngleArray, inv_P0_extent2);
    H_tilda = p * H_f;
    
    % Obtain the prediction for the nonlinear implicit measurement
    % model, h(x,y) = 0
    iMeasModelPredicted = P2 * R_from_G_to_L * (iMeas3D_G - predCenter)...
        - meanS * H_tilda * predExtent_P2;
    predMeas_P2(1+(i-1)*2 : i*2) = iMeasModelPredicted;
    
    % Obtain the covariance of the current measurement
    % Definition: iCovMeas = k(uk,uk), uk: argument of the current measurement
    iCovMeas = compute_GP_covariance(angleMeas_L, angleMeas_L, paramGP{2});
    % Definition: iR_f = k(uk,uk) - k(uk, uf) * inv(K(uf,uf)) * k(uf, uk)
    iR_f = iCovMeas - covMeasBasis * inv_P0_extent2 * covMeasBasis';
    % Compute the total covariance of this specific measurement
    iR =  stdMeasGP^2 * eye(2) + stdS^2*p*iR_f*p' + stdS^2*H_tilda...
        *(predExtent_P2*predExtent_P2')* H_tilda';
    % Plug this matrix into the overall covariance matrix
    measCov_P2(((i-1)*2+1):i*2, ((i-1)*2+1):i*2) = iR;
    
    %% Obtain the partial derivatives of the measurement model
    [H_center, H_velocity, H_quatDev, H_angVel, H_extent, H_meas] = compute_linearized_measurement_matrix...
        (iMeas3D_G, predState, prevEstQuat, P2, predExtent_P2...
        , paramGP{2}, basisAngleArray, inv_P0_extent2, meanS, eps);
    
    dH_dX_P2(1+(i-1)*2:i*2, :) = [H_center H_velocity H_quatDev H_angVel...
        zeros(2, numBasisAngles) H_extent zeros(2, numBasisAngles)];
    
    dH_dY_P2(1+(i-1)*2:i*2, :) = H_meas;
end

% Plane 3
for i = 1:curNumMeas
    iMeas3D_G = transpose(measArray(i, :));
    
    % Compute the following variables to exploit in measurement model
    [p, H_f, covMeasBasis, angleMeas_L] = compute_measurement_model_components...
        (iMeas3D_G, predState, prevEstQuat, P3, paramGP{3}...
        , basisAngleArray, inv_P0_extent3);
    H_tilda = p * H_f;
    
    % Obtain the prediction for the nonlinear implicit measurement
    % model, h(x,y) = 0
    iMeasModelPredicted = P3 * R_from_G_to_L * (iMeas3D_G - predCenter)...
        - meanS * H_tilda * predExtent_P3;
    predMeas_P3(1+(i-1)*2 : i*2) = iMeasModelPredicted;
    
    % Obtain the covariance of the current measurement
    % Definition: iCovMeas = k(uk,uk), uk: argument of the current measurement
    iCovMeas = compute_GP_covariance(angleMeas_L, angleMeas_L, paramGP{3});
    % Definition: iR_f = k(uk,uk) - k(uk, uf) * inv(K(uf,uf)) * k(uf, uk)
    iR_f = iCovMeas - covMeasBasis * inv_P0_extent3 * covMeasBasis';
    iR = stdMeasGP^2 * eye(2) + stdS^2*p*iR_f*p' + stdS^2*H_tilda...
        *(predExtent_P3*predExtent_P3')* H_tilda';
    % Plug this matrix into the overall covariance matrix
    measCov_P3(((i-1)*2+1):i*2, ((i-1)*2+1):i*2) = iR;
    
    %% Obtain the partial derivatives of the measurement model
    [H_center, H_velocity, H_quatDev, H_angVel, H_extent, H_meas] = compute_linearized_measurement_matrix...
        (iMeas3D_G, predState, prevEstQuat, P3, predExtent_P3...
        , paramGP{3}, basisAngleArray, inv_P0_extent3, meanS, eps);
    
    dH_dX_P3(1+(i-1)*2:i*2, :) = [H_center H_velocity H_quatDev H_angVel...
        zeros(2, numBasisAngles) zeros(2, numBasisAngles) H_extent];
    
    dH_dY_P3(1+(i-1)*2:i*2, :) = H_meas;
end


% Put the measurements in a column the form: [x1; y1;...; xn; yn]
predMeas =  [predMeas_P1; predMeas_P2; predMeas_P3];
% Measurement Update
dH_dX = [dH_dX_P1; dH_dX_P2; dH_dX_P3];
dH_dY = [dH_dY_P1; dH_dY_P2; dH_dY_P3];
measCov = blkdiag(measCov_P1, measCov_P2, measCov_P3);
kalmanGain = predStateCov * dH_dX' / (dH_dX*predStateCov*dH_dX' + measCov);
estState = predState + kalmanGain * (0 - predMeas);
estStateCov = (eye(numStates) - kalmanGain*dH_dX) * predStateCov;
estStateCov = (estStateCov + estStateCov')/2; % To make it symmetric (eliminating numeric errors)

% Quaternion Treatment
estQuatDev = estState(7:9);
estQuat = apply_quat_deviation(estQuatDev, prevEstQuat);
estState(7:9) = zeros(3,1);     % Reset the orientation deviation
end

function [ H_center, H_velocity, H_quatDev, H_angVel, H_extent, H_meas] = compute_linearized_measurement_matrix...
    (meas3D_G, state, quatPrev, projectionMatrix, extent, kernelParameters, basisAngleArray, inv_P0_extent, meanS, eps)

center = state(1:3);
quatDev = state(7:9);
quat = apply_quat_deviation(quatDev, quatPrev);

R_from_G_to_L = rotation_matrix_from_global_to_local(quat);

%% Obtain the measurement prediction by the original nonlinear implicit meas model
[p, H_f, ~, ~] = compute_measurement_model_components(meas3D_G, state, quatPrev...
    , projectionMatrix, kernelParameters, basisAngleArray, inv_P0_extent);
measModelOriginal = projectionMatrix * R_from_G_to_L * (meas3D_G - center)...
    - meanS * p * H_f * extent;

%% Calculate partial derivative wrt target extent
H_extent = - meanS * p * H_f;     % This partial derivative is taken analytically.

%% Calculate partial derivative wrt linear positions
H_center = zeros(2,3);
for i = 1:3
    stateIncremented = state;   % Initialize the vector
    stateIncremented(i) = stateIncremented(i) + eps;    % increment the state
    centerIncremented = stateIncremented(1:3);          % incremented center
    
    % Compute the output of original measurement model for the incremented state vector
    [p, H_f, ~, ~] = compute_measurement_model_components(meas3D_G, stateIncremented, quatPrev...
        , projectionMatrix, kernelParameters, basisAngleArray, inv_P0_extent);
    measModelForIncState= projectionMatrix * R_from_G_to_L * (meas3D_G - centerIncremented)...
        - meanS * p * H_f * extent;
    
    % Compute numeric derivative
    H_center(:, i) = (measModelForIncState - measModelOriginal) / eps;
end

%% Calculate partial derivative wrt linear velocities
H_velocity = zeros(2,3); % Measurement model does not depend on linear velocities

%% Calculate partial derivative wrt quaternion deviations
H_quatDev = zeros(2,3);
for i = 7:9
    stateIncremented = state;   % Initialize the vector
    stateIncremented(i) = stateIncremented(i) + eps;
    quatIncremented = apply_quat_deviation(stateIncremented(7:9), quatPrev);
    
    % Compute the output of original measurement model for the incremented state vector
    [p, H_f, ~, ~] = compute_measurement_model_components(meas3D_G, stateIncremented, quatPrev...
        , projectionMatrix, kernelParameters, basisAngleArray, inv_P0_extent);
    
    % Compute the rotation matrix for the incremented quaternion
    R_from_G_to_L_Incremented = rotation_matrix_from_global_to_local(quatIncremented);
    
    % Compute the output of original measurement model for the incremented state vector
    measModelForIncState = projectionMatrix * R_from_G_to_L_Incremented * (meas3D_G - center)...
        - meanS * p * H_f * extent;
    
    % Compute numeric derivative
    H_quatDev(:,i-6) = (measModelForIncState - measModelOriginal) / eps;
end

%% Calculate partial derivative wrt angular velocities
H_angVel = zeros(2,3); % Measurement model does not depend on linear velocities

%% Calculate partial derivatives wrt measurement
H_meas = projectionMatrix * R_from_G_to_L;

end

function [diffUnitVector, H_f, covMeasBasis, angle] = compute_measurement_model_components...
    (meas3D_G, state, quatPrev, projectionMatrix, paramGP, basisAngleArray, inv_P0_extent)

center = state(1:3);
quatDev = state(7:9);
quat = apply_quat_deviation(quatDev, quatPrev); % Apply predicted orientation deviation

% Compute the following variables to exploit in measurement model
R_from_G_to_L = rotation_matrix_from_global_to_local(quat);

% Transform the measurements to the local frame
meas3D_L = R_from_G_to_L * (meas3D_G - center);

% Project the local measurement onto the corresponding projection plane
measProjected = projectionMatrix * meas3D_L;

diffVector = measProjected;
diffVectorMag = norm(diffVector);
if diffVectorMag == 0
    diffVectorMag = eps;
end
diffUnitVector = diffVector / diffVectorMag;

angle = atan2(diffVector(2), diffVector(1));  % angle of the projected measurement on the projection plane
angle = mod(angle, 2*pi);

covMeasBasis = compute_GP_covariance(angle, basisAngleArray, paramGP);

H_f = covMeasBasis * inv_P0_extent; % GP model relating the extent and the radial
end

function [inv_P0_extent] = compute_inverse_GP_covariance(basisAngleArray, paramGP, eps)
P0_extent = compute_GP_covariance(basisAngleArray, basisAngleArray, paramGP);
P0_extent = P0_extent + eps*eye(size(basisAngleArray,1));       % to prevent numerical errors thrown in matrix inversion
chol_P0_extent = chol(P0_extent);
inv_chol_P0_extent = inv(chol_P0_extent);
inv_P0_extent = inv_chol_P0_extent * inv_chol_P0_extent';
end

function [F, Q] = compute_rotational_dynamics(w, std, T)
% Inputs:
%              w:       Angular velocity
%              std:     Standard deviation of the angular velocity
%              T:       Sampling time

% Dummy variables
wNorm = norm(w);
if wNorm == 0
    wNorm = 1e-3;
end

S = skew_symmetric_matrix(-w);
c = cos(1/2*T*wNorm);
s = sin(1/2*T*wNorm);
expMat = eye(3,3) + s/wNorm*S + (1-c)/wNorm^2*S^2;

% Construct the state transition matrix
F = [expMat  T*expMat ; zeros(3,3)  eye(3,3)];

% Construct the process noise covariance matrix
G11 = T*eye(3,3) + 2/wNorm^2*(1-c)*S + 1/wNorm^2*(T-2/wNorm*s)*S^2;
G12 = 1/2*T^2*eye(3,3) + 1/wNorm^2*(4/wNorm*s-2*T*c)*S + 1/wNorm^2*(1/2*T^2+2/wNorm*T*s+4/wNorm^2*(c-1))*S^2;
G21 = zeros(3,3);
G22 = T * eye(3,3);
G = [G11 G12; G21 G22];

B = G * [zeros(3,3); eye(3,3)];

% cov = std^2 * diag([0 0 1]);
cov = std^2 * eye(3,3);
Q = B * cov * B';
end

function [out] = skew_symmetric_matrix(x)
out = [0  -x(3)  x(2);...
    x(3)  0  -x(1);...
    -x(2)  x(1)  0];
end

function [qOut] = apply_quat_deviation(a, qIn)
qa = 1/ sqrt(4+norm(a)) * [a; 2];   % It depends on the Rodrigues parametrization
qOut = quat_product(qa, qIn);

qOut = qOut/ norm(qOut);    % Normalize due to numeric inprecisions
end

function [out] = quat_product(q1, q2)
q1V = q1(1:3);
q1S = q1(4);
q2V = q2(1:3);
q2S = q2(4);

out = [q1S*q2V + q2S*q1V - cross(q1V,q2V);...
    q1S*q2S - q1V'*q2V];
end