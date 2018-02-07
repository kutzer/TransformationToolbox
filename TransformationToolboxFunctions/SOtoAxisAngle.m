function [Axis,Angle] = SOtoAxisAngle(R)
% SOtoAxisAngle converts a 2D or 3D rotation matrix (element of SO(2) or 
%   SO(3) to axis/angle.
%
%   Note: This method includes an extension to n-dimensional rotations
%
%   [Axis,Angle] = SOtoAxisAngle(R)
%
%   M. Kutzer 22Jan2016, USNA

% Updates
%   02Feb2016 - Updated to include n-dimensional axis/angle.
%   07Feb2018 - Updated to replace SO check with a warning.

% TODO - update documentation to reflect n-dimensional axis/angle
% TODO - address negative eigenvalue issues of logm for larger than 3x3

%% Check inputs
narginchk(1,1);
[bin,msg] = isSO(R);
if ~bin
    warning('SOtoAxisAngle:NotSO',...
        ['Input must be a valid 2D or 3D rotation matrix.\n',...
        ' -> %s'],msg);
end

%% Calculate axis angle
N = size(R,1);
switch N
    case 2
        Angle = atan2(R(2),R(1));
        Axis = sign(Angle);
        Angle = abs(Angle);
    case 3
        r = vrrotmat2vec(R);
        Axis = r(1:3);
        Angle = r(4);
    otherwise
        %error('SOtoAxisAngle:NotSO',...
        %    'Input must be a valid 2D or 3D rotation matrix.');
        
        % TODO - address negative eigenvalue issues of logm
        r = logm(R);
        v = vee(r,'fast');
        Angle = norm(v);
        Axis = transpose( v./Angle );
end