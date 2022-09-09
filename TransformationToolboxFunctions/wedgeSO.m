function M = wedgeSO(v,varargin)
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
%       M - MxM real skew-symmetric matrix element of so(M)
%
%   See also wedge, vee, veeSO, soBasis, isSkewSymmetric, veeSE, wedgeSE.
%
%   M. Kutzer 04Jan2017, USNA

% Update(s)
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast [UNUSED]

%% Check inputs
narginchk(1,3);

%{
% UNUSED
% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms
%}

%% Calculate wedge
M = wedge(v,varargin(:));