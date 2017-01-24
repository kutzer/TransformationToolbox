function v = vee(M,options)
% VEE converts an NxN skew-symmetric matrix into a vector defined by the
% basis elements of so(N) (the Lie algebra associated with SO(n), sometimes 
% referred to as "little so"). 
%   v = VEE(M) converts an NxN skew-symmetric matrix "M" into a Mx1 vector 
%   "v".
%       2 x 2 matrix -> $v \in \mathbb{R}^1$ (2D rotations)
%       3 x 3 matrix -> $v \in \mathbb{R}^3$ (3D rotations)
%       4 x 4 matrix -> $v \in \mathbb{R}^6$
%       5 x 5 matrix -> $v \in \mathbb{R}^10$
%       ...
%       N x N matrix -> $v \in \mathbb{R}^M$
%
%   v = VEE(___,'fast') converts without checking the NxN matrix for 
%   skew-symmetry. 
%
%   See also veeSO, wedge, wedgeSO, soBasis, isSkewSymmetric, veeSE, 
%            wedgeSE.
%
%   M. Kutzer 10Oct2014, USNA

% Updates
%   03Feb2016 - Documentation update and added check for isSkewSymmetric.m
%   04Jan2017 - Updated documentation
%   24Jan2017 - Replaced skew-symmetric error 
%% Default options
narginchk(1,2);
if nargin < 2
    options = '';
end

%% Check M
switch lower(options)
    case 'fast'
        % Do not check for skew-symmetry
    otherwise
        if exist('isSkewSymmetric','file')
            [bin,msg] = isSkewSymmetric(M);
            if ~bin
                warning('"M" must be skew-symmetric.\n\t -> %s',msg);
            end
        else
            warning('The function "isSkewSymmetric" was not found. Skipping skew-symmetry check.');
        end
end

%% Calculate v
if size(M,1) == 3
    v(1,1) = M(3,2);
    v(2,1) = M(1,3);
    v(3,1) = M(2,1);
    return
end

if size(M,1) == 2
    v(1,1) = M(2,1);
    return
end

%% $M \in so(N)$
n = size(M,1);
e = soBasis(n);
m = numel(e);
for idx = 1:m
    [i,j] = find(e{idx} == 1);
    v(idx,1) = M(i,j);
end

end