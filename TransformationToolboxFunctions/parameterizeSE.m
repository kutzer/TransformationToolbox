function [t_a2b,r_a2b] = parameterizeSE(H_a2b,method)
% PARAMETERIZESE parameterizes one or more elements of SE(N) using a
% user-specified method
%   [t_a2b,r_a2b] = parameterizeSE(H_a2b)
%
%   [t_a2b,r_a2b] = parameterizeSE(H_a2b,method)
%
%   Input(s)
%       H_a2b  - M-element cell array containing elements of SE(N)
%       method - [OPTIONAL] character or string argument defining method
%                for pareterization
%
%   Output(s)
%       t_a2b - KxM array containing translation elements of
%               parameterization
%       r_a2b - PxM array containing rotation elements of parameterization
%
%   M. Kutzer, 02Nov2022, USNA

%% Check input(s)
narginchk(1,2)

% TODO - add ZERO, true/false, etc matching TransformationToolbox tools

if nargin < 2
    method = 'decoupled';
end

% TODO - check method

% Account for non-cell input
if ~iscell(H_a2b)
    H_a2b = {H_a2b};
end

%% Define size parameters
M = numel(H_a2b);
n = size(H_a2b{1},1);

%% Define method
switch lower(method)
    case 'coupled'
        func = @(H)coupled(H,n);
    case 'decoupled'
        func = @(H)decoupled(H,n);
    otherwise
        error('Method is not currently supported.')
end

%% Define additional parameters
[t_a2b,r_a2b] = func(H_a2b{1});
K = size(t_a2b,1);
P = size(r_a2b,1);

%% Parameterize
t_a2b = zeros(K,M);
r_a2b = zeros(P,M);
for i = 1:M
    [t_a2b(:,i),r_a2b(:,i)] = func(H_a2b{i});
end

end

%% Internal functions
function [t,r] = decoupled(H,n)
localFast = true;
t = H(1:(n-1),n);
R = H(1:(n-1),1:(n-1));
r = veeSO( logSO(R,localFast),localFast );
end

function [t,r] = coupled(H,n)
localFast = true;
v = veeSE( logSE(H,localFast),localFast );
k = numel(v);
r = v(1:(k-(n-1)));
t = v((k-(n-1)+1):end);
end