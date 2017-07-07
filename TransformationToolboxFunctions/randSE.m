function H = randSE(r,dim)
% randSE creates a random homogeneous transformation element of SE(n) using
% pseudorandom numbers.
%   randSE creates a random 3D homogeneous transformation with a
%   translation norm bounded by [0,1]
%   randSE(r) creates a random 3D homogeneous transformation with a
%   translation norm bounded by [0,r]
%   randSE(r,dim) creates a random element of SE(n) with a translation norm
%   bounded by [0,r]
%
%   See also randSO, wedge, vee, wedgeSO, veeSO
%
%   M. Kutzer 07July2017, USNA

%Updates
% 

%% Check inputs
if nargin < 1
    r = 1;
end

if nargin < 2
    dim = 3;
end

if dim < 2
    error('Transformation must be of at least dimension 2.');
end

%% Define random transformation
R = randSO(dim);
T_0 = 2*rand(dim,1) - repmat(1,dim,1);
T_hat = T_0./norm(T_0);
T = r*rand(1,1)*T_hat;

H = R;
H(:,end+1) = T;
H(end+1,end) = 1;