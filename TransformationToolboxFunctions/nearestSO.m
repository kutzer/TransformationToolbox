function R = nearestSO(M)
% NEARESTSO finds the nearest element of the special orthogonal group to a
% given matrix.
%   R = NEARESTSO(M) calculates the nearest rotation matrix R to a given
%   matrix M such that the square of the Frobenius norm is minimized.
%    -> R is the argument that minimizes \| M - R \|_{F}^{2}
%
%   Input(s)
%       M - an nxn matrix (n must match a valid dimension for SO(N))
%   
%   Output(s)
%       R - nxn element of SO(N)
%
%   References:
%       "Finding the Nearest Orthonormal Matrix," accessed 25Sep2020.
%       http://people.csail.mit.edu/bkph/articles/Nearest_Orthonormal_Matrix.pdf
%
%   M. Kutzer, 25Sep2020, USNA 

%% Check Input(s)
narginchk(1,1);

% TODO - check dimension of M

%% Calculate R
% TODO - check if (transpose(M)*M)^(-1/2) is possible
R = M*(transpose(M)*M)^(-1/2);