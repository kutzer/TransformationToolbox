function H = SOtoSE(R)
% SOtoSE converts an element of SO(n) to SE(n).
%   H = SOtoSE(R) 
%
%   M. Kutzer, 28Jun2021, USNA

%% Check input(s)
narginchk(1,1);
[bin,msg] = isSO(R);
if ~bin
    warning('Specified rotation matrix is not valid: %s',msg);
end

%% Append row/column
n = size(R,1);
H = eye(n+1);
H(1:n,1:n) = R;
