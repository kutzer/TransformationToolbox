function R = expSO(r,varargin)
% EXPSO calculates the matrix exponential of an element of the special
% orthogonal group.
%   R = EXPSO(r) calculates the matrix exponential using the general
%   formulation for SO(2), and Rodrigues's formula for SO(3). SO(N > 3)
%   uses the matrix exponential. 
%
%   R = EXPSO(___,ZERO)
%
%   R = EXPSO(___,fast)
%
%   Input(s)
%       r      - NxN element of so(N)
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
%
%   M. Kutzer, 26Jan2022, USNA

% Update(s)
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast

%% Default options
ZERO = [];
fast = false;

%% Check inputs
narginchk(1,3);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

if ~fast
    % TODO - check h for valid so(N) properties
end

% TODO - Check bad inputs (cellOut)

%% Convert to vector
k = vee(r,ZERO,fast);

%% Calculate axis/angle
Angle = norm(k);
if Angle == 0
    R = eye(size(r,1));
else
    Axis = transpose( k./Angle );
    R = AxisAngletoSO(Axis,Angle);
end