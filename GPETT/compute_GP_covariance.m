function [ covMatrix ] = compute_GP_covariance(argArray1, argArray2, paramGP)
% This function computes the GP covariance according to the specified parameters.
% Output:
%           covMatrix:      The covariance matrix computed by the GP kernel.
% Inputs:
%           argArray1:       Each row is a pair of spherical angles, i.e., = [azimuth elevation]
%           argArray2:       Each row is a pair of spherical angles, i.e., = [azimuth elevation]
%           paramGP:         Parameters of the underlying GP

% Extract parameters
kernelType = paramGP{1};
stdPrior = paramGP{2};
stdRadius  = paramGP{3};
scaleLength = paramGP{4};
isObjectSymmetric = paramGP{5};

switch kernelType
    
    case 1
        % Periodic squared exponential kernel               
        diffMatrix = compute_diffence_matrix(argArray1, argArray2);
        
        if isObjectSymmetric
            % Periodic with pi (it imposes the symmetry of the object)
            covMatrix = stdPrior^2 * exp(-sin(diffMatrix).^2 / (2*scaleLength^2)) + stdRadius^2;
        else
            % Covariance periodic with 2pi
            covMatrix = stdPrior^2 * exp(-diffMatrix.^2 / (2*scaleLength^2)) + stdRadius^2;
        end
        
    case 2
        % Matern 3_2 kernel
        diffMatrix = compute_diffence_matrix(argArray1, argArray2);        
        covMatrix = stdPrior^2 * (1+sqrt(3)/scaleLength*diffMatrix) .* exp(-sqrt(3)/scaleLength*diffMatrix);
        
    case 3
        % Matern 5_2 kernel
        diffMatrix = compute_diffence_matrix(argArray1, argArray2);        
        covMatrix = stdPrior^2 * (1+sqrt(5)/scaleLength*diffMatrix...
            + 5/(3*scaleLength^2)*diffMatrix.^2) .* exp(-sqrt(5)/scaleLength*diffMatrix);
end
end


function [ diffMatrix ] = compute_diffence_matrix( inp1, inp2 )
% This function calculates the difference matrix regarding two inputs.
% inp1, inp2: Each row is of the form [azimuthAngle elevationAngle]

% Produce a grid structure to be able to compute each angle between
% elements of two inputs.
len1 = size(inp1, 1); % The number of spherical angles, : [azi_1 elv_1;... azi_N elv_N]
len2 = size(inp2, 1); % The number of spherical angles: : [azi_1 elv_1;... azi_N elv_N]

% Extract azimuth and elevation angles from the inputs
aziArray1 = inp1(:,1);
elvArray1 = inp1(:,2);
aziArray2 = inp2(:,1);
eleArray2 = inp2(:,2);

% Create matrices from these vectors to form grid structure
aziGrid1 = repmat(aziArray1, [1, len2]);
elvGrid1 = repmat(elvArray1, [1, len2]);
aziGrid2 = repmat(transpose(aziArray2), [len1, 1]);
elvGrid2 = repmat(transpose(eleArray2), [len1, 1]);

% Put the matrices into vectors to feed to the angle computation procedure
aziArray1 = aziGrid1(:);
elvArray1 = elvGrid1(:);
aziArray2 = aziGrid2(:);
eleArray2 = elvGrid2(:);

% Prepare the arrays for the following function
sphrCoorArray1 = [aziArray1 elvArray1];
sphrCoorArray2 = [aziArray2 eleArray2];

% Calculate the the difference (spherical angle)
diffArray = compute_angle_btw_sphr_coor(sphrCoorArray1, sphrCoorArray2);

% Transform difference array back in matrix form
diffMatrix = reshape(diffArray, len1, len2);
end

function [ angle ] = compute_angle_btw_sphr_coor( sphrCoor1, sphrCoor2 )
% The function computes the angle between two Spherical coordinates.
% sphrCoor1, sphrCoor2 are in the form of : [azimuthAngle elevationAngle]. Note that
% azimuthAngle and elevation angle can also be column vectors.

azi1 = sphrCoor1(:, 1);
elv1 = sphrCoor1(:, 2);
azi2 = sphrCoor2(:, 1);
elv2 = sphrCoor2(:, 2);

% Calculate trigonometric multiplications
multCosAzi = cos(azi1) .* cos(azi2);
multSinAzi = sin(azi1) .* sin(azi2);
multSinElv = sin(elv1) .* sin(elv2);
multCosElv = cos(elv1) .* cos(elv2);

angle = real(acos(multCosElv.*multCosAzi + multCosElv.*multSinAzi + multSinElv));
end
