function M = wedge(v)
% WEDGE converts an Nx1 vector into a skew-symmetric matrix 
% for all N \in {1,3,6,10,15,21,28,...}.
%   M = WEDGE(v) calculates a skew-symmetric matrix from elements of v. 
%   Note that the number of elements in v must correspond to the non-zero
%   upper-triangular elements of a real, skew-symmetric matrix:
%       2 x 2 matrix -> $v \in \mathbb{R}^1$ (2D rotations)
%       3 x 3 matrix -> $v \in \mathbb{R}^3$ (3D rotations)
%       4 x 4 matrix -> $v \in \mathbb{R}^6$
%       5 x 5 matrix -> $v \in \mathbb{R}^10$
%       etc.
%
%   See also wedgeSO, vee, veeSO, soBasis, isSkewSymmetric, veeSE, wedgeSE.
%
%   M. Kutzer 09Oct2014, USNA

%Updates
%   27Feb2015 - Updated to include N-dimensional wedge (or hat) operator.
%   03Feb2016 - Updated documentation
%   04Jan2017 - Updated documentation

%% Check inputs
narginchk(1,1);
m = numel(v);
v = reshape(v,m,[]);

%% $M \in so(3)$
if m == 3        
    M = [    0, -v(3),  v(2);...
          v(3),     0, -v(1);...
         -v(2),  v(1),     0];
    return
end

%% $M \in so(2)$
if m == 1
    M = [   0, -v(1);...
         v(1),     0];
    return
end

%% $M \in so(N)$
%m = numel(v);
%m = n*(n-1)/2;
n = (8*m + 1)^(1/2)/2 + 1/2;
if n ~= round(n)
    error('wedge:BadVector',...
        ['Invalid dimension for input vector.\n',...
        ' -> The specified input vector must correspond to the number\n',...
        '    of non-zero upper-triangular elements of a real,\n',...
        '    skew-symmetric matrix.'])
end
e = soBasis(n);
M = zeros(size(e{1}));
for i = 1:m
    M = M + v(i)*e{i};
end 

end