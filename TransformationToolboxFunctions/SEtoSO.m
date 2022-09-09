function [R,d] = SEtoSO(H,varargin)
% SEtoSO isolates the element of SO(N) from SE(N)
%   R = SEtoSO(H)
%
%   [R,d] = SEtoSO(H)
%
%   ... = SEtoSO(___,ZERO)
%
%   ... = SEtoSO(___,fast)
%
%   Input(s)
%       H      - (N+1)x(N+1) element of SE(N)
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
%       R - NxN element of SO(N)
%       d - Nx1 array containing translation associated with SE(N)
%
%   M. Kutzer, 27Jan2022, USNA

% Update(s)
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast

%% Set defaults
ZERO = [];
fast = false;

%% Check input(s)
narginchk(1,3);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms

if ~fast
    [bin,msg] = isSE(H,ZERO);
    if ~bin
        warning('Specified rigid body transform matrix is not valid: %s',msg);
    end
end

%% Parse outputs
N = size(H,1);
R = H(1:(N-1),1:(N-1));
d = H(1:(N-1),end);