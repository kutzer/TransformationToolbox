function [bin,msg] = isSE(M)
% ISSE checks a matrix to see if it is an element of the special Euclidean
% group
%   ISSE(M) checks an nxn matrix "M" for the properties of the special 
%   Euclidean group. If M is an element of the special Euclidean
%   group, this function returns "1", "0" is returned otherwise.
%
%   msg - message describing property that is violated
%   
%   See also hgtransform triad showTriad hideTriad isSO
%
%   M. Kutzer, 13May2015, USNA

%% Check dimensions
d = size(M);
if numel(d) ~= 2 || (d(1) ~= d(2))
    msg = 'Matrix is not NxN.';
    bin = 0;
    return
end

%% Check if matrix is real
if ~isreal(M)
    msg = 'Matrix is not real.';
    bin = 0;
    return
end

%% Check rotation matrix
n = size(M,1);
[bin,msg] = isSO(M(1:n-1,1:n-1));
if bin == 0
    msg = sprintf('Rotation %s',msg);
    return
end

%% Otherwise
bin = 1;
msg = [];