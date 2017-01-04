function [bin,msg] = isSO(M)
% ISSO checks a matrix to see if it is an element of the special orthogonal
% group
%   ISSO(M) checks an nxn matrix "M" for the properties of the special 
%   orthogonal group. If M is an element of the special orthogonal
%   group, this function returns "1", "0" is returned otherwise.
%
%   msg - message describing property that is violated
%   
%   See also isSE
%
%   M. Kutzer 12May2015, USNA

% Updates
%   03Jan2017 - Updated to relax constraints on the determinant, mutual
%               orthogonality, and unit length rows and columns.

%% Check dimensions
d = size(M);
if numel(d) ~= 2 || (d(1) ~= d(2))
    msg = 'Matrix is not NxN.';
    bin = 0;
    return
end
n = d(1);

%% Check if matrix is real
if ~isreal(M)
    msg = 'Matrix is not real.';
    bin = 0;
    return
end

%% Check for determinant of 1
ZERO = 2e3*eps(class(M));
detM = det(M);
if ~isZero(detM-1,ZERO)
    msg = sprintf('Matrix has a determinant of %.15f.',detM);
    bin = 0;
    return
end
    
%% Check for orthogonality/inverse property
ZERO = 2e3*eps(class(M));
I = M*M';
if ~isZero(I-eye(size(I)),ZERO)
    msg = sprintf('Matrix has columns/rows that are not mutually orthogonal.\n');
    bin = 0;
    for i = 1:size(I,1)
        for j = 1:size(I,2)
            msg = [msg,sprintf('\tI(%d,%d) = %.15f\n',i,j,I(i,j))];
        end
    end
    return
end

%% Check unit vector length of columns/rows
ZERO = 2e3*eps(class(M));
magM = sqrt(sum(M.^2,1));
if ~isZero(magM-ones(size(magM)),ZERO)
    msg = sprintf('Matrix has columns/rows that are not unit length.\n');
    for i = 1:numel(magM)
        msg = [msg,sprintf('\t|M(:,%d)| = %.15f\n',i,magM(i))];
    end
    bin = 0;
    return
end

%% Otherwise
msg = [];
bin = 1;

