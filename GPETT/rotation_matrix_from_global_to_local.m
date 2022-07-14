function [ rotMatrix ] = rotation_matrix_from_global_to_local( quat )
% The function computes the rotation matrix from the global frame to the Local (body)
% frame.
% Note that the rotation matrix is to multiply the vector from LEFT,
% i.e., xL  =  rotMatrix  *  xG.

quat = quat(:);     % Makes sure the input is a column matrix

q = quat(1:3);
q4 = quat(4);

% Rotation matrix from global to local R_from_G_to_L
rotMatrix = transpose((q4^2-q'*q)*eye(3,3) + 2*(q*q') - 2*q4*skew_symmetric_matrix(q));
end

function [out] = skew_symmetric_matrix(x)
out = [0  -x(3)  x(2);...
    x(3)  0  -x(1);...
    -x(2)  x(1)  0];
end
