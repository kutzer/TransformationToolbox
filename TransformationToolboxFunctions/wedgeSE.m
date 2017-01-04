function M = wedgeSE(v)
% WEDGESE converts an Mx1 vector into an element of se(n) (the Lie algebra 
% associated with SE(n), sometimes referred to as "little se") for all 
% M \in {3,6,10,15,21,28,...}.
%   M = WEDGESE(v) calculates an element of se(n) from elements of v.
%   Note that the number of elements in v must correspond to the number of
%   basis elements of se(n).
%       3 x 3 matrix -> $v \in \mathbb{R}^3$ (2D rigid body motion)
%       4 x 4 matrix -> $v \in \mathbb{R}^6$ (3D rigid body motion)
%       5 x 5 matrix -> $v \in \mathbb{R}^10$
%       6 x 6 matrix -> $v \in \mathbb{R}^15$
%       ...
%
%   See also veeSE, seBasis, veeSO, wedgeSO.
%
%   M. Kutzer, 04Jan2017, USNA

%% Check inputs
narginchk(1,1);
p = numel(v);
v = reshape(v,p,[]);

%% $M \in so(N)$
m = p - (8*p + 1)^(1/2)/2 + 1/2;
n = (8*m + 1)^(1/2)/2 + 1/2;
if n ~= round(n)
    error('wedge:BadVector',...
        ['Invalid dimension for input vector.\n',...
        ' -> The specified input vector must correspond to the number\n',...
        '    of basis elements of an associated se(n).'])
end
e = seBasis(n);
M = zeros(size(e{1}));
for i = 1:p
    M = M + v(i)*e{i};
end 