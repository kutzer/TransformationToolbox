function d = distanceSO(R1,R2,varargin)
% DISTANCESO calculates the distance between elements of SO(N) using the 
% "\phi_6 method" [1] to quantify distance between elements of SO(N).
%       d = ||log(R1 R2^{\top})|| or
%       d = |vee( log(R1 R2^{\top})|
%
%   d = DISTANCESO(R1,R2) 
%
%   d = DISTANCESO(___,ZERO)
%
%   d = DISTANCESO(___,FAST)
%
%   Input(s)
%       R1      - R1 is an NxN matrix element of SO(N)
%       R2      - R2 is an NxN matrix element of SO(N)
%       ZERO   - [OPTIONAL] positive scalar value that is
%                sufficiently close to zero to be assumed zero
%                (e.g. ZERO = 1e-8). If a "ZERO" is not specified,
%                a default of ZERO = 1e-8 is used.
%       fast   - [OPTIONAL] true/false logical value indicating
%                whether to skip checking a specified property or
%                properties. If "fast" is not specified, a default of
%                fast = false is used.
%           fast = true    - Skip checking property or properties
%           fast = [false] - Check property or properties
%
%   Output(s)
%       d - scalar value quantifying distance between R1 and R2
%           d \in [0,\pi]
%
% Reference(s)
% [1] Huynh, Du Q. "Metrics for 3D rotations: Comparison and analysis." 
%     Journal of Mathematical Imaging and Vision 35.2 (2009): 155-164.
%
%   M. Kutzer, 08Sep2021, USNA

% Update(s)
%   09Nov2022 - Updated to use parseVarargin_ZERO_fast and updated
%               documentation

% TODO - add alternative metrics, see reference.
% TODO - address warning issue?

%% Default options
ZERO = 1e-8;
fast = false;

%% Check inputs
narginchk(1,4);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

if ~fast
    if ~isSO(R1,ZERO)
        error('R1 must be an NxN matrix element of SO(N).');
    end
    if ~isSO(R2,ZERO)
        error('R2 must be an NxN matrix element of SO(N).');
    end
end

N1 = size(R1);
N2 = size(R2);
if nnz(N1 == N2) ~= 2
    error('R1 and R2 must be NxN matrix elements of SO(N).');
end

%% Calculate 
fastLocal = true;
d = norm( logSO(R1*R2.',fastLocal),2 );
