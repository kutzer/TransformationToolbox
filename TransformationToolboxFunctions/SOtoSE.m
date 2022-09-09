function H = SOtoSE(R,varargin)
% SOtoSE converts an element of SO(n) to SE(n).
%   H = SOtoSE(R) 
%
%   H = SOtoSE(___,ZERO)
%
%   H = SOtoSE(___,fast)
%
%   Input(s)
%       R - NxN element of SO(N)
%       ZERO   - [OPTIONAL] positive scalar value that is
%                sufficiently close to zero to be assumed zero
%                (e.g. ZERO = 1e-8). If a "ZERO" is not specified,
%                a default of ZERO = [] is used.
%       fast   - [OPTIONAL] true/false logical value indicating
%                whether to skip checking a specified property or
%                properties. If "fast" is not specified, a default of
%                fast = false is used.
%           fast = true    - Skip checking property or properties
%           fast = [false] - Check property or properties
%
%   Output(s)
%       H - (N+1)x(N+1) element of SE(N)
%
%   M. Kutzer, 28Jun2021, USNA

% Update(s)
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast and documentation
%               updates

% TODO - consider adding d (translation) as an optional input to match
%        SEtoSO

%% Set defaults
ZERO = [];
fast = false;

%% Check input(s)
narginchk(1,3);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms

if ~fast
    [bin,msg] = isSO(R,ZERO);
    if ~bin
        warning('Specified rotation matrix is not valid: %s',msg);
    end
end

%% Append row/column
n = size(R,1);
H = eye(n+1);
H(1:n,1:n) = R;
