function [bin,msg] = isSE(M,ZERO)
% ISSE checks a matrix to see if it is an element of the special Euclidean
% group
%   ISSE(M) checks an nxn matrix "M" for the properties of the special 
%   Euclidean group. If M is an element of the special Euclidean
%   group, this function returns "1", "0" is returned otherwise.
%
%   ISSE(M,ZERO) specifies a value sufficiently close to zero that
%   properties are checked against
%
%   [bin,msg] = ISSE(___) returns both a binary (t/f) value and, if 
%   applicable, a string describing why the matrix is not an element of 
%   SE(N)
%
%   Input(s)
%       M    - NxN numeric array
%       ZERO - positive value that is sufficiently close to zero to be
%              assumed zero (e.g. ZERO = 1e-8). If ZERO is not specified or
%              if ZERO = [], a default value is used.
%
%   Output(s)
%       bin - logical value describing whether the matrix is an element of 
%             SE(N).
%       msg - message describing property that is violated. This parameter
%             is [] if bin is "true".
%
%   See also hgtransform triad showTriad hideTriad isSO
%
%   M. Kutzer, 13May2015, USNA

% Update(s)
%   18Nov2021 - Updated to check last row for [0,...,0,1]
%   18Nov2021 - Added optional ZERO input
%   18Nov2021 - Updated documentation
%   15Sep2021 - Account for NxN where N < 3
%   02Apr2025 - Updated to make "bin" always logical

%% Check input(s)
narginchk(1,2);
if nargin < 2
    ZERO = [];
end

%% Check dimensions
d = size(M);
if numel(d) ~= 2 || (d(1) ~= d(2))
    msg = 'Matrix is not NxN.';
    bin = false;
    return
end
n = d(1);

if n < 3
    msg = 'Matrix must be at least 3x3.';
    bin = false;
    return
end

%% Check if matrix is real
if ~isreal(M)
    msg = 'Matrix is not real.';
    bin = false;
    return
end

%% Check last row for [0,...,0,1]
v0 = M(end,1:n-1);
v1 = M(end,end);
if ~isZero(v0,ZERO) || ~isZero(v1-1,ZERO)
    msg = 'Matrix does not contain a last row of [';
    for i = 1:n-1
        msg = sprintf('%s0, ',msg);
    end
    msg = sprintf('%s1].',msg);
    bin = false;
    return
end


%% Check rotation matrix
[bin,msg] = isSO(M(1:n-1,1:n-1),ZERO);
if ~bin
    msg = sprintf('Rotation %s',msg);
    return
end

%% Otherwise
bin = true;
msg = [];