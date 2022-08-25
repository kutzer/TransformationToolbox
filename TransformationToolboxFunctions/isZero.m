function bin = isZero(M,ZERO)
%ISZERO checks each element of an array to see if it zero. If all elements
%are zero, isZero returns a "1" and "0" otherwise.
%   ISZERO(M) checks each element of an arbitrary array M against a
%   value for zero defined using 10*eps(class(M))
%
%   ISZERO(M,ZERO) checks each element of an arbitrary array M against a
%   specified value for zero. If no value for zero is specified, a default
%   value is specified using the class of M and the spacing of floating
%   point numbers for that associated class (using eps.m).
%
%   Note: NaN values are treated as non-zero
%
%   bin = ISZERO(___)
%
%   Input(s)
%       M    - numeric array of arbitrary dimensions
%       ZERO - positive value that is sufficiently close to zero to be
%              assumed zero (e.g. ZERO = 1e-8). If ZERO is not specified or
%              if ZERO = [], ZERO is defined as 10*eps(class(M))
%
%   Output(s)
%       bin - value describing whether all values contained in M are within
%             +/- ZERO of 0
%
%   See also zeroFPError eps
%
%   M. Kutzer 13May2015, USNA

% Updates
%   22Jan2016 - Updated to speed up checking matrices
%   18Nov2021 - Added optional ZERO input
%   18Nov2021 - Updated documentation
%   25Aug2022 - Updated to account for NaN values

%% Set default zero
defaultZERO = 10*eps(class(M));
if nargin < 2
    ZERO = defaultZERO;
end
if isempty(ZERO)
    ZERO = defaultZERO;
end

%% Check if a matrix is effectively zero or NaN
% Check for values that are effectively zero
BIN = abs(M) > ZERO;
% Check for nan values
BIN = BIN | isnan(M);

if nnz(BIN) > 0
    bin = 0;
else
    bin = 1;
end