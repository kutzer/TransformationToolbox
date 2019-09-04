function invH = invSE(H)
% INVSE Calculates the inverse of an element of the Special Euclidean group
% using the properties of rotation matrices.
%   
%   See also invSO
%
%   M. Kutzer, 15May2015, USNA

% Updates:
%   04Sep2019 - Added details to error message
%   04Sep2019 - Replaced error with warning message

%% Check input
[bin,msg] = isSE(H);
if ~bin
    warning('invSE:NotSE','Input must be a valid member of the Special Euclidean group.\n\t-> %s',msg);
end

%% Calculate inverse
n = size(H,1);
R = H(1:n-1, 1:n-1);
V = H(1:n-1,n);

invH = eye(n);
invH(1:n-1,1:n-1) = transpose(R);
invH(1:n-1,n) = -transpose(R)*V;