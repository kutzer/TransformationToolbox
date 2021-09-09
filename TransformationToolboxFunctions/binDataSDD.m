function [idx,nBins] = binDataSDD(x,n)
% BINDATASDD clusters data into bins with bin bounds defined using the
% standard deviation of the difference found using the sorted data set.
%   idx = binDataSDD(x) uses a default value of n = 2
%   idx = binDataSDD(x,n)
%
%   Input(s)
%       x - N element array containing real values
%       n - scalar, positive, real value defining multiples of the standard 
%           deviation of difference to define the bin boundaries.
%               68.27%, n = 1
%               95.45%, n = 2
%               99.73%, n = 3
%               ...
%
%   Output(s)
%       idx   - index values associated with bins
%           x_1 = x(idx == 1); % return values in the 1st bin (smallest values)
%           ...
%           x_n = x(idx == nBins); % return values in last bin (largest values)
%       nBins - total number of bins found
%
%   M. Kutzer, 09Sep2021, USNA

%% Parse input(s)
sizex = size(x);
x = reshape(x,1,numel(x));

if ~isreal(x)
    error('x must contain only real values.');
end

if nargin < 2
    n = 2;
end

if numel(n) > 1
    error('n must be a scalar value.');
end

if ~isreal(n)
    error('n must be a real value.');
end

if n < 0
    error('n must be a positive value.');
end

%% Cluster data
[sx,Ix] = sort(x);
dsx = diff(sx);
sigma = std(dsx);
idx = find(dsx > n*sigma);

vals = [sx(idx); sx(idx+1)];
dmu = mean(vals,1);

edges = [0,dmu,max(x)];
idx = discretize(x,edges);

%% Package output(s)
idx = reshape(idx,sizex);
nBins = max(idx);