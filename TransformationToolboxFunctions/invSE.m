function invH = invSE(H)
% INVSE Calculates the inverse of an element of the Special Euclidean group
% using the properties of rotation matrices.
%   
%   See also invSO
%
%   (c) M. Kutzer 15May2015, USNA

%% Check input
if ~isSE(H)
    error('Input must be a valid member of the Special Euclidean group.');
end

%% Calculate inverse
n = size(H,1);
R = H(1:n-1, 1:n-1);
V = H(1:n-1,n);

invH = eye(n);
invH(1:n-1,1:n-1) = transpose(R);
invH(1:n-1,n) = -transpose(R)*V;