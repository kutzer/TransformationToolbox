function invR = invSO(R,varargin)
% INVSO Calculates the inverse of an element of the Special Orthogonal 
% group using the properties of rotation matrices (i.e. the inverse is the 
% transpose).
%   invR = INVSO(R)
%
%   invR = INVSO(___,ZERO)
%
%   invR = INVSO(___,fast)
%
%   Input(s)
%       R    - NxN array element of SO(N)
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not   
%              specified, a default value of ZERO = [] is used.
%       fast - [OPTIONAL] true/false logical value indicating whether to
%              skip checking SE(N). Choosing fast = true ignores specified
%              ZERO. 
%                fast = true    - Skip checking if H \in SE(N)
%                fast = [false] - Check if H \in SE(N)
%
%   Output(s)
%       invR - NxN element of SO(N) representing the matrix inverse of R
%
%   See also invSE
%
%   M. Kutzer, 03Feb2016, USNA

% Updates:
%   04Sep2019 - Added details to error message
%   04Sep2019 - Replaced error with warning message
%   06Sep2022 - Updated to include ZERO and fast optional inputs
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast

%% Default options
ZERO = [];
fast = false;

%% Check inputs
narginchk(1,3);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms

%% Check R
if ~fast
    [bin,msg] = isSO(R,ZERO);
    if ~bin
        warning('invSO:NotSO','Input must be a valid member of the Special Orthogonal group.\n\t-> %s',msg);
    end
end

%% Calculate inverse
invR = transpose( R );