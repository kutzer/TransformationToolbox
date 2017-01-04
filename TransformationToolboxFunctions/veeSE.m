function v = veeSE(M,options)
% VEESE converts an NxN element of se(n) (the Lie algebra associated with
% SE(n), sometimes referred to as "little se") into a vector defined by the
% basis elements of se(N).
%   v = VEESE(M) converts an NxN matrix "M" into an Mx1 vector "v".
%       3 x 3 matrix -> $v \in \mathbb{R}^3$ (2D rigid body motion)
%       4 x 4 matrix -> $v \in \mathbb{R}^6$ (3D rigid body motion)
%       5 x 5 matrix -> $v \in \mathbb{R}^10$
%       6 x 6 matrix -> $v \in \mathbb{R}^15$
%       ...
%       N x N matrix -> $v \in \mathbb{R}^M$
%
%   v = VEESE(___,'fast') converts without checking the (N-1)x(N-1) matrix
%   for skew-symmetry.
%
%   See also wedgeSE, seBasis, veeSO, wedgeSO.
%
%   M. Kutzer, 04Jan2017, USNA

%% Default options
narginchk(1,2);
if nargin < 2
    options = '';
end

%% Check M
N = size(M,1);
switch lower(options)
    case 'fast'
        % Do not check (N-1)x(N-1) for skew-symmetry
    otherwise
        if exist('isSkewSymmetric','file')
            if ~isSkewSymmetric(M(1:(N-1),1:(N-1)))
                error('"M" must be skew-symmetric.');
            end
        else
            warning('The function "isSkewSymmetric" was not found. Skipping skew-symmetry check.');
        end
end

%% $M \in se(N)$
n = N-1;
e = seBasis(n);
m = numel(e);
for idx = 1:m
    [i,j] = find(e{idx} == 1);
    v(idx,1) = M(i,j);
end