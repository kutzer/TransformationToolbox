function v = nCross(varargin)
% NCROSS calculates the n-dimensional cross product following the the
% definition of the wedge product (exterior algebra).
%   v = nCross(v1,v2,v3,...,vN)
%
%   Input(s)
%       v1 - M-element array
%       v2 - M-element array
%       v3 - M-element array
%       ...
%       vN - M-element array
%
%       NOTE: number of arrays supplied must be M-1 (one less vector than
%             the vector's dimension).
%
%   Output(s)
%       v - M-element array
%   
%   M. Kutzer 12Mar2015, USNA

% Updates
%   24Feb2022 - Updated to correct 2D cross product, added input(s) check, 
%               and updated documentation

% TODO - add references


%% Define number of vectors and dimensions
N = numel(varargin);     % Total number of vectors supplied
nD = numel(varargin{1}); % Dimension of vectors

%% Check inputs
% Check for appropriate number of vectors
if N ~= nD-1
    error('For an %dD vector cross product, %d vectors must be supplied.',nD,nD-1);
end

% Check vector dimensions
for i = 2:N
    if numel(varargin{i}) ~= nD
        error('All vectors must be the same dimension.');
    end
end

%% Combine into matrix form
E = zeros(nD);
for i = 2:(N+1)
    E(i,:) = reshape(varargin{i-1},1,nD);
end

%% Select signs
if nD == 2 
    e = [1, -1];
end 

if nD == 3
    e = [1, -1, 1];
end

if nD == 4
    e = [1,-1,1,-1];
end

%TODO - confirm the correct signs for e(i) for nD > 4
if nD > 4
    for i = 1:nD
        e(i) = (-1)^(i-1);
    end
end

%% Calculate cross product
E(1,:) = e;
for i = 1:nD
    M = [E(2:end,1:(i-1)),E(2:end,(i+1):end)];
    v(i,1) = e(i)*det(M);
end

%% Output size matching the first input vector
v = reshape(v,size(varargin{1}));