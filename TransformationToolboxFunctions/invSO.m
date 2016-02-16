function invR = invSO(R)
% INVSO Calculates the inverse of an element of the Special Orthogonal 
% group using the properties of rotation matrices (i.e. the inverse is the 
% transpose).
%   
%   See also invSE
%
%   M. Kutzer 03Feb2016, USNA


%% Check input
if ~isS0(R)
    error('Input must be a valid member of the Special Orthogonal group.');
end

invR = transpose( R );