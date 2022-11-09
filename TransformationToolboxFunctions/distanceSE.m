function [dist_t,dist_r] = distanceSE(H_a2o,H_b2o,varargin)
% DISTANCESE calculates the between elements of SE(N). Rotation distances
% are calculated using the "\phi_6 method" [1] to quantify distance between
% elements of SO(N). Distance in translation is calculated using the norm.
%
%   [dist_t,dist_r] = DISTANCESE(H_a2o,H_b2o) 
%
%   [dist_t,dist_r] = DISTANCESE(___,ZERO)
%
%   [dist_t,dist_r] = DISTANCESE(___,FAST)
%
%   Input(s)
%       H_a2o  - H_b2o is an (N+1)x(N+1) matrix element of SE(N)
%       H_b2o  - H_b2o is an (N+1)x(N+1) matrix element of SE(N)
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
%       dist_t - scalar value quantifying translational distance between 
%                H_a2o and H_b2o
%       dist_r - scalar value quantifying rotational distance between H_a2o
%                and H_b2o. Note that dist_r \in [0,\pi]
%
% Reference(s)
% [1] Huynh, Du Q. "Metrics for 3D rotations: Comparison and analysis." 
%     Journal of Mathematical Imaging and Vision 35.2 (2009): 155-164.
%
%   M. Kutzer, 09Nov2022, USNA

%% Default options
ZERO = 1e-8;
fast = false;

%% Check inputs
narginchk(1,4);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

if ~fast
    if ~isSE(H_a2o,ZERO)
        error('H_a2o must be an (N+1)x(N+1) matrix element of SE(N).');
    end
    if ~isSE(H_b2o,ZERO)
        error('H_b2o must be an (N+1)x(N+1) matrix element of SE(N).');
    end
end

N1 = size(H_a2o);
N2 = size(H_b2o);
if nnz(N1 == N2) ~= 2
    error('H_a2o and H_b2o must be (N+1)x(N+1) matrix element of SE(N).');
end

%% Calculate relative rotation and translation
fastLocal = true;
H_a2b = H_a2o*invSE(H_b2o,fastLocal);

N = size(H_a2b,1) - 1;
R_a2b = H_a2b(1:N,1:N);
d_a2b = H_a2b(1:N,N+1);

%% Calculate translational distance
dist_t = norm(d_a2b);

%% Calculate rotational distance [1]
dist_r = norm( logSO(R_a2b,fastLocal),2 );