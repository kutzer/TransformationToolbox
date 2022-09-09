function r = logSO(R,varargin)
% LOGSO calculates the matrix natural log of an element of the special 
% orthogonal group.
%   r = LOGSO(R) calculates the matrix natural log using the general 
%   formulation for SO(2), and Rodrigues's formula for SO(3). SO(N > 3)
%   uses the matrix natural log.
%
%   r = LOGSO(___,ZERO)
%
%   r = LOGSO(___,fast)
%
%   Input(s)
%       R - NxN element of SO(N)
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not   
%              specified, a default value is used.
%       fast - [OPTIONAL] true/false logical value indicating whether to
%              skip checking SO(N). Choosing fast = true ignores specified
%              ZERO. 
%                fast = true    - Skip checking if R \in SO(N)
%                fast = [false] - Check if R \in SO(N)
%
%   Output(s)
%       r - NxN element of so(N)
%
%   See also expSO SOtoAxisAngle AxisAngletoSO logSE expSE
%
%   M. Kutzer 08Jan2016, USNA

% Updates
%   26Jan2022 - Replaced logm usage entirely
%   26Jan2022 - Updated documentation
%   06Sep2022 - Updated to include ZERO and fast optional inputs
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast

%% Set default values
ZERO = [];
fast = false;

%% Check input(s)
narginchk(1,3);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms

%% Calculate logSO
[Axis,Angle] = SOtoAxisAngle(R,ZERO,fast);
v = Axis*Angle;
r = wedge(v);
