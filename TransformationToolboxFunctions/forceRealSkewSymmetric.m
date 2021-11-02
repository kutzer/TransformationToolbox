function B = forceRealSkewSymmetric(A)
% FORCEREALSKEWSYMMETRIC forces a matrix to be real skew-symmetric by
% zeroing the diagonal and averaging complimentary off-diagonal terms.
%   B = FORCEREALSKEWSYMMETRIC(A)
%
%   Input(s)
%       A - NxN real matrix
%
%   Output(s)
%       B - NxN skew-symmetric matrix
%
%   M. Kutzer, 02Nov2021, USNA

%% Check input(s)
[m,n,o] = size(A);
if m ~= n
    error('Input matrix must be square.');
end
if o ~= 1
    error('Input matrix must be NxN.');
end

%% Keep real parts of A
if ~isreal(A)
    warning('Matrix has imaginary parts');
    A = real(A);
end

%% Define skew-symmetric 
uA = triu(A, 1); % Upper triangular
lA = tril(A,-1); % Lower triangular

uB = (uA - lA.')./2;

B = uB - uB.';