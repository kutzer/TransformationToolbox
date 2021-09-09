function d = distanceSO(R1,R2,varargin)
% DISTANCESO calculates the between elements of SO(N).
%   d = DISTANCESO(R1,R2) uses the \phi_6 method to quantify distance
%   between elements of SO(N)
%       d = ||log(R1 R2^{\top})|| or
%       d = |vee( log(R1 R2^{\top})|
%
%   Input(s)
%       R1 - R1 is an NxN matrix element of SO(N)
%       R2 - R2 is an NxN matrix element of SO(N)
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

% TODO - add alternative metrics, see reference.
% TODO - address warning issue?

%% Check input(s)
N1 = size(R1);
N2 = size(R2);
if numel(N1) ~= 2
    error('R1 must be an NxN matrix element of SO(N).');
end
if numel(N2) ~= 2
    error('R2 must be an NxN matrix element of SO(N).');
end
if nnz(N1 == N2) ~= 2
    error('R1 and R2 must be NxN matrix elements of SO(N).');
end

N = N1(1);

%% Calculate 
warning off
d = norm( logSO(R1*R2.'),2 );
warning on