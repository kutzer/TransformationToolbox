function bin = isZero(M,ZERO)
%ISZERO checks each element of an array to see if it zero. If all elements
%are zero, isZero returns a "1" and "0" otherwise.
%   ISZERO(M,ZERO) checks each element of an arbitrary array M against a
%   specified value for zero. If no value for zero is specified, a default
%   value is specified using the class of M and the spacing of floating
%   point numbers for that associated class (using eps.m).
%
%   See also zeroFPError eps
%
%   M. Kutzer 13May2015, USNA

% Updates
%   22Jan2016 - Updated to speed up checking matrices

%% Set default zero
if nargin < 2
    ZERO = 10*eps(class(M));
end

%% Check if a matrix is effectively zero
BIN = abs(M) > ZERO;
if sum(reshape(BIN,1,[])) > 0
    bin = 0;
else
    bin = 1;
end