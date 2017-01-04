function M = wedgeSO(v)
% WEDGESO converts an Nx1 vector into a skew-symmetric matrix 
% for all N \in {1,3,6,10,15,21,28,...}.
%   M = WEDGESO(v) calculates a skew-symmetric matrix from elements of v. 
%   Note that the number of elements in v must correspond to the non-zero
%   upper-triangular elements of a real, skew-symmetric matrix:
%       2 x 2 matrix -> $v \in \mathbb{R}^1$ (2D rotations)
%       3 x 3 matrix -> $v \in \mathbb{R}^3$ (3D rotations)
%       4 x 4 matrix -> $v \in \mathbb{R}^6$
%       5 x 5 matrix -> $v \in \mathbb{R}^10$
%       etc.
%
%   See also wedge, vee, veeSO, soBasis, isSkewSymmetric, veeSE, wedgeSE.
%
%   M. Kutzer 04Jan2017, USNA

%% Check inputs
narginchk(1,1);

%% Calculate wedge
M = wedge(v);