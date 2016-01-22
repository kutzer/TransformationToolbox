function [bin,msg] = isSO(M)
% isSO checks a matrix to see if it is an element of the special orthogonal
% group
%   isSO(M) checks an nxn matrix "M" for the properties of the special 
%   orthogonal group. If M is an element of the special orthogonal
%   group, this function returns "1", "0" is returned otherwise.
%
%   msg - message describing property that is violated
%   
%   See also isSE
%
%   (c) M. Kutzer 12May2015, USNA

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

%% Check for determinant of 1
detM = det(M);
if ~isZero(detM-1)
    msg = sprintf('Matrix has a determinant of %.15f.',detM);
    bin = 0;
    return
end
    
%% Check for orthogonality/inverse property
I = M*M';
if ~isZero(I-eye(size(I)))
    msg = 'Matrix has columns/rows that are not mutually orthogonal.';
    bin = 0;
    return
end

%% Check unit vector length of columns/rows
magM = sqrt(sum(M.^2,1));
if ~isZero(magM-ones(size(magM)))
    msg = 'Matrix has columns/rows that are not unit length.';
    bin = 0;
    return
end

%% Otherwise
msg = [];
bin = 1;

