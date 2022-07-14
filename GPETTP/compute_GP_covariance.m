function [ covMatrix ] = compute_GP_covariance(argArray1, argArray2, paramGP)
% This function computes the GP covariance according to the specified parameters.
% Output:
%           covMatrix:      The covariance matrix computed by the GP kernel.
% Inputs:
%           argArray1:       Each row is a polar angle.
%           argArray2:       Each row is a polar angle.
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
            % Covariance periodic with pi (it imposes the symmetry assumpation of the object)
            covMatrix = stdPrior^2 * exp(-sin(diffMatrix).^2 / (2*scaleLength^2)) + stdRadius^2;
        else
            % Covariance periodic with 2*pi
            covMatrix = stdPrior^2 * exp(-2 * sin(diffMatrix/2).^2 / scaleLength^2) + stdRadius^2;
        end                
end
end

function [ diffMatrix ] = compute_diffence_matrix( inp1, inp2 )
% This function calculates the angle matrix regarding two inputs.

% Produce a grid structure to be able to compute each angle between
% elements of two inputs.
len1 = size(inp1, 1); % The number of polar angles in the first input
len2 = size(inp2, 1); % The number of polar angles in the second input

% Create matrices from these vectors to form grid structure
grid1 = repmat(inp1, [1, len2]);
grid2 = repmat(transpose(inp2), [len1, 1]);

% Compute the difference matrix to be used for covariance matrix calculation
diffMatrix = grid1-grid2;
end
