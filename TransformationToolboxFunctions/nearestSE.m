function H = nearestSE(M)
% NEARESTSE finds the nearest element to the special Euclidean group given
% a matrix.
%   H = nearestSE(M) 
%
%   Input(s)
%       M - NxN array containing a matrix element "close" to being an 
%           element of SE(N-1). N must be greater than 2.
%   
%   Output(s)
%       H - NxN array describing the "closest" element of SE(N-1)
%
%   Notes:
%       (1) This method does not properly map reflections to SE(3).
%       (2) This method does not properly map matrices containing
%           "rotation" (upper left 3x3 block) elements with a determinant 
%           of 0. 
%
%   Example(s):
%       % Define random value of SE(3)
%       H_tru = randSE;
%
%       % Alter value to an estimate "close" to SE(3)
%       H_est = H_tru;
%       H_est(1:3,1:3) = 0.8*H_est(1:3,1:3);
%
%       % Recover nearest element of SE(3)
%       H_rec = nearestSE(H_est);
%       
%       % Compare results
%       H_dif = H_tru*invSE(H_rec)
%
%   Failed Example(s)
%       % (1) Reflection Example
%       % Define reflection
%       H_ref = eye(4);
%       H_ref(3,3) = -1;
%       H_rec = nearestSE(H_ref)
%       [tf,msg] = isSE(H_rec)
%
%       % (2) Zero Determinant Example
%       H_ref = eye(4);
%       H_ref(3,3) = 0;
%       H_rec = nearestSE(H_ref)
%       [tf,msg] = isSE(H_rec)
%
%   M. Kutzer, 25Sep2020, USNA

% Updates
%   02Sep2022 - Updated documentation to include examples and limitations

%% Check input(s)
narginchk(1,1);

[m,n] = size(M);
if ~ismatrix(M) || m ~= n || m < 3
    error(['Input matrix must be square matrix larger than 2x2. ',...
        'Current matrix input is %dx%d.'],m,n);
end

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

% TODO - consider checking final result