function v = veeSO(M,varargin)
% VEESO converts an NxN skew-symmetric matrix into a vector defined by the
% basis elements of so(N) (the Lie algebra associated with SO(n), sometimes 
% referred to as "little so"). 
%   v = VEESO(M) converts an NxN skew-symmetric matrix "M" into a Mx1  
%   vector "v".
%       2 x 2 matrix -> $v \in \mathbb{R}^1$ (2D rotations)
%       3 x 3 matrix -> $v \in \mathbb{R}^3$ (3D rotations)
%       4 x 4 matrix -> $v \in \mathbb{R}^6$
%       5 x 5 matrix -> $v \in \mathbb{R}^10$
%       ...
%       N x N matrix -> $v \in \mathbb{R}^M$
%
%   v = VEESO(___,'fast') converts without checking the NxN matrix for 
%   skew-symmetry. 
%
%   See also vee, wedge, wedgeSO, soBasis, isSkewSymmetric, veeSE, wedgeSE.
%
%   M. Kutzer 04Jan2017, USNA


%% Calculate vee
v = vee(M,varargin(:));