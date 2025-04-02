function [bin,msg] = isSO(M,ZERO)
% ISSO checks a matrix to see if it is an element of the special orthogonal
% group
%   ISSO(M) checks an nxn matrix "M" for the properties of the special 
%   orthogonal group. If M is an element of the special orthogonal
%   group, this function returns "1", "0" is returned otherwise.
%
%   ISSO(M,ZERO) specifies a value sufficiently close to zero that
%   properties are checked against
%
%   [bin,msg] = ISSO(___) returns both a binary (t/f) value and, if 
%   applicable, a string describing why the matrix is not an element of 
%   SO(N)
%
%   Input(s)
%       M    - NxN numeric array
%       ZERO - positive value that is sufficiently close to zero to be
%              assumed zero (e.g. ZERO = 1e-8). If ZERO is not specified or
%              if ZERO = [], a default value is used.
%
%   Output(s)
%       bin - logical value describing whether the matrix is an element of 
%             SO(N).
%       msg - message describing property that is violated. This parameter
%             is [] if bin is "true".
%   
%   See also isSE
%
%   M. Kutzer, 12May2015, USNA

% Updates
%   03Jan2017 - Updated to relax constraints on the determinant, mutual
%               orthogonality, and unit length rows and columns
%   07Feb2018 - Updated to actively calculate ZERO based on test condition
%   18Nov2021 - Added optional ZERO input
%   18Nov2021 - Updated documentation
%   02Apr2025 - Updated to make "bin" always logical

%% Check input(s)
narginchk(1,2);
if nargin < 2
    ZERO = [];
end

%% Define ZERO scaling term
ZERO_scale = 1e1;

%% Check dimensions
d = size(M);
if numel(d) ~= 2 || (d(1) ~= d(2))
    msg = 'Matrix is not NxN.';
    bin = false;
    return
end
%n = d(1);

%% Check if matrix is real
if ~isreal(M)
    msg = 'Matrix is not real.';
    bin = false;
    return
end

%% Check for determinant of 1
detM = det(M);

% Set ZERO
%defaultZERO = 2e3*eps(class(M));
defaultZERO = ZERO_scale * max([eps(detM),1e-6]);
if ~isempty(ZERO)
    ZERO_i = ZERO;
else
    ZERO_i = defaultZERO;
end

if ~isZero(detM-1,ZERO_i)
    msg = sprintf('Matrix has a determinant of %.15f.',detM);
    bin = false;
    return
end
    
%% Check for orthogonality/inverse property
I = M*M';

% Set ZERO
%defaultZERO = 2e3*eps(class(M));
defaultZERO = ZERO_scale * max([ max(reshape(eps(I),1,[])), max(reshape(eps(eye(size(I))),1,[])) ]);
if ~isempty(ZERO)
    ZERO_i = ZERO;
else
    ZERO_i = defaultZERO;
end

if ~isZero(I-eye(size(I)),ZERO_i)
    msg = sprintf('Matrix has columns/rows that are not mutually orthogonal.\n');
    msg = [msg,sprintf('\tConsider updating ZERO from %e to %e\n',ZERO_i,max(abs(reshape(I-eye(size(I)),1,[]))))];
    bin = false;
    for i = 1:size(I,1)
        for j = 1:size(I,2)
            msg = [msg,sprintf('\t\tI(%d,%d) = %.15f\n',i,j,I(i,j))];
        end
    end
    return
end

%% Check unit vector length of columns/rows
magM = sqrt(sum(M.^2,1));

% Set ZERO
%defaultZERO = 2e3*eps(class(M));
defaultZERO = ZERO_scale * max([ max(eps(magM)), max(eps(ones(size(magM)))) ]);
if ~isempty(ZERO)
    ZERO_i = ZERO;
else
    ZERO_i = defaultZERO;
end

if ~isZero(magM-ones(size(magM)),ZERO_i)
    msg = sprintf('Matrix has columns/rows that are not unit length.\n');
    for i = 1:numel(magM)
        msg = [msg,sprintf('\t|M(:,%d)| = %.15f\n',i,magM(i))];
    end
    bin = false;
    return
end

%% Otherwise
msg = [];
bin = true;
