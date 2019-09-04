function invR = invSO(R)
% INVSO Calculates the inverse of an element of the Special Orthogonal 
% group using the properties of rotation matrices (i.e. the inverse is the 
% transpose).
%   
%   See also invSE
%
%   M. Kutzer, 03Feb2016, USNA

% Updates:
%   04Sep2019 - Added details to error message
%   04Sep2019 - Replaced error with warning message

%% Check input
[bin,msg] = isSO(R);
if ~bin
    warning('invSO:NotSO','Input must be a valid member of the Special Orthogonal group.\n\t-> %s',msg);
end

invR = transpose( R );