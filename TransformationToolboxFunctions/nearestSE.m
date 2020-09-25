function H = nearestSE(M)
% NEARESTSE finds the nearest element to the special Euclidean group given
% a matrix.
%   H = nearestSE(M) 
%
%   M. Kutzer, 25Sep2020, USNA

% TODO - finish documentation

%% Check input(s)
narginchk(1,1);

% TODO - check dimension of M

%% Break up matrix
n = size(M,1);
% Rotation Element
M_R = M(1:(n-1),1:(n-1));
% Translation Element
M_d = M((1:n-1),n);

%% Define rotation
R = nearestSO(M_R);

%% Define translation
% TODO - confirm that this is best practice.
d = real(M_d);

%% Build transformation
H = eye(n,n);
H(1:(n-1),1:(n-1)) = R;
H(1:(n-1),n) = d;

