function invH = invSE(H,varargin)
% INVSE Calculates the inverse of an element of the Special Euclidean group
% using the properties of rotation matrices.
%   invH = INVSE(H)
%
%   invH = LOGSE(___,ZERO)
%
%   invH = LOGSE(___,fast)
%
%   Input(s)
%       H    - (N+1)x(N+1) array element of SE(N)
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not   
%              specified, a default value is used.
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

%% Default options
ZERO = [];
fast = false;

%% Check inputs
narginchk(1,3);

if nargin > 1
    for i = 1:numel(varargin)
        switch lower( class(varargin{i} ))
            case 'logical'
                fast = varargin{i};
            otherwise
                if numel(varargin{i}) ~= 1
                    error('Numeric optional inputs must be scalar values.');
                end

                if varargin{i} == 0 || varargin{i} == 1
                    fast = logical(varargin{i});
                else
                    ZERO = varargin{i};
                end
        end
    end
end

% Check zero
if ZERO < 0
    error('ZERO value must be greater or equal to 0.')
end

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