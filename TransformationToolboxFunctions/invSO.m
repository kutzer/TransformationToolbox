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
%              specified, a default value is used.
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

%% Check R
if ~fast
    [bin,msg] = isSO(R);
    if ~bin
        warning('invSO:NotSO','Input must be a valid member of the Special Orthogonal group.\n\t-> %s',msg);
    end
end

%% Calculate inverse
invR = transpose( R );