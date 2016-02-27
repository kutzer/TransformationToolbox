function v = nCross(varargin)
% NCROSS calculates the n-dimensional cross product following the the
% definition of the wedge product (exterior algebra).
%
%   M. Kutzer 12Mar2015, USNA

% TODO - check number of inputs
% TODO - update documentation

%% Define vectors for cross product
N = numel(varargin);
n = numel(varargin{1});

E = zeros(n);
for i = 2:(N+1)
    E(i,:) = reshape(varargin{i-1},n,1);
end

%% select signs
if n == 2 
    e = 1;
end 

if n == 3
    e = [1, -1, 1];
end

if n == 4
    e = [1,-1,1,-1];
end
%TODO - actually figure out the correct signs for e!
if n > 4
    for i = 1:n
        e(i) = (-1)^(i-1);
    end
end

%% calculate cross product
E(1,:) = e;
for i = 1:n
    M = [E(2:end,1:(i-1)),E(2:end,(i+1):end)];
    v(i,1) = e(i)*det(M);
end
    