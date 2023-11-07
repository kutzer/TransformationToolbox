function H = randSE(r,dim)
% randSE creates a random, n-dimensional homogeneous rigid body 
% transformation element of SE(n).
%
%   Use Options:
%
%   OPTION 1: 
%   H_a2b = randSE creates a random 3D homogeneous transformation with a
%   translation norm bounded by [0,1] and rotation is generated using
%   randSO.
%
%   H_a2b = randSE(r) creates a random 3D homogeneous transformation with a
%   translation norm bounded by [0,r] and rotation is generated using
%   randSO.
%
%   H_a2b = randSE(r,dim) creates a random N-D homogeneous transformation 
%   (N = dim) with a translation norm bounded by [0,r] and rotation is 
%   generated using randSO.
%
%   OPTION 2 (TODO - Complete Option 2)
%   H_a2b = randSE(H_a2b_cov) generates a random sample from a 
%   multivariate normal distribution (see mvnrnd) assuming a zero mean and 
%   using the specified rigid body transformation covariance (see covSE).
%
%   H_a2b = randSE(H_a2b_mu,H_a2b_cov) generates a random sample from a 
%   multivariate normal distribution (see mvnrnd) given a specified mean 
%   (see meanSE) and covariance (see covSE).
%
%   H_a2b = randSE(___,n) allows the user to specify the number of samples
%   returned. Samples will be returned in a cell array. 
%
%   [H_a2b,x_a2b,k_a2b] = randSE(___)
%
%   See also randSO, wedge, vee, wedgeSO, veeSO, wedgeSE, veeSE, meanSE,
%   covSE
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