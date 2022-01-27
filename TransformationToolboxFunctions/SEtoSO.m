function [R,d] = SEtoSO(H)
% SEtoSO isolates the element of SO(N) from SE(N)
%   R = SEtoSO(H)
%   [R,d] = SEtoSO(H)
%
%   Input(s)
%       H - (N+1)x(N+1) element of SE(N)
%
%   Output(s)
%       R - NxN element of SO(N)
%       d - Nx1 array containing translation associated with SE(N)
%
%   M. Kutzer, 27Jan2022, USNA

%% Check input(s)
narginchk(1,1);
[bin,msg] = isSE(H);
if ~bin
    warning('Specified rigid body transform matrix is not valid: %s',msg);
end

%% Parse outputs
N = size(H,1);
R = H(1:(N-1),1:(N-1));
d = H(1:(N-1),end);