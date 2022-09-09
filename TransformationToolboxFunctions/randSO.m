function R = randSO(dim)
% randSO creates a random rotation element of SO(n) using pseudorandom
% numbers.
%   randSO creates a random 3D rotation matrix
%   randSO(n) creates a random element of SO(n)
%
%   See also randSE, wedge, vee, wedgeSO, veeSO
%
%   M. Kutzer 07July2017, USNA

% Updates
%   

%% Check inputs
if nargin < 1
    dim = 3;
end

if dim < 2
    error('Rotation must be at least 2x2.');
end

%% Define vector dimension
n = dim;
m = ((2*(n - 1/2))^2 - 1)/8;

%% Create random rotation vector
v = 2*pi*rand(m,1);
V = wedgeSO(v);
R = expm(V);