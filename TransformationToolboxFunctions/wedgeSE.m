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
%   Input(s)
%       v - N-element array (see notes above regarding valid values of N)
%
%   Input(s) [Unused, added for standard input syntax]
%       ZERO   - [OPTIONAL] positive scalar value that is
%                sufficiently close to zero to be assumed zero
%                (e.g. ZERO = 1e-8). If a "ZERO" is not specified,
%                a default of ZERO = [] is used.
%       fast   - [OPTIONAL] true/false logical value indicating
%                whether to skip checking a specified property or
%                properties. If "fast" is not specified, a default of
%                fast = false is used.
%
%   Output(s)
%       M - MxM real matrix element of se(M)
%
%   See also veeSE, seBasis, veeSO, wedgeSO.
%
%   M. Kutzer, 04Jan2017, USNA

% Update(s)
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast [UNUSED]

%% Check inputs
narginchk(1,3);
p = numel(v);
v = reshape(v,p,[]);

%{
% UNUSED
% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms
%}

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