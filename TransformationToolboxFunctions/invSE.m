function invH = invSE(H,varargin)
% INVSE Calculates the inverse of an element of the Special Euclidean group
% using the properties of rotation matrices.
%   invH = INVSE(H)
%
%   invH = INVSE(___,ZERO)
%
%   invH = INVSE(___,fast)
%
%   Input(s)
%       H    - (N+1)x(N+1) array element of SE(N)
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not   
%              specified, a default of ZERO = [] is used.
%       fast - [OPTIONAL] true/false logical value indicating whether to
%              skip checking SE(N). Choosing fast = true ignores specified
%              ZERO. 
%                fast = true    - Skip checking if H \in SE(N)
%                fast = [false] - Check if H \in SE(N)
%
%   Output(s)
%       invH - (N+1)x(N+1) element of SE(N) representing the matrix inverse
%              of H
%
%   See also invSO
%
%   M. Kutzer, 15May2015, USNA

% Updates:
%   04Sep2019 - Added details to error message
%   04Sep2019 - Replaced error with warning message

% Updates
%   02Nov2021 - Updated to include "fast"
%   08Sep2022 - Updated to include ZERO 
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast

%% Default options
ZERO = [];
fast = false;

%% Check inputs
narginchk(1,3);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms

%% Check H
if ~fast
    [bin,msg] = isSE(H,ZERO);
    if ~bin
        warning('invSE:NotSE','Input must be a valid member of the Special Euclidean group.\n\t-> %s',msg);
    end
end

%% Calculate inverse
n = size(H,1);
R = H(1:n-1, 1:n-1);
V = H(1:n-1,n);

% Initialize inverse 
invH = eye(n);
% Inverse rotation
invH(1:n-1,1:n-1) = transpose(R);
% Inverse translation
invH(1:n-1,n) = -transpose(R)*V;