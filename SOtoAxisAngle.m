function [Axis,Angle] = SOtoAxisAngle(R)
% SOtoAxisAngle converts a 2D or 3D rotation matrix (element of SO(2) or 
%   SO(3) to axis/angle.
%
%   [Axis,Angle] = SOtoAxisAngle(R)
%
%   M. Kutzer 22Jan2016, USNA

%% Check inputs
narginchk(1,1);
[bin,msg] = isSO(R);
if ~bin
    error('SOtoAxisAngle:NotSO',...
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
        error('SOtoAxisAngle:NotSO',...
            'Input must be a valid 2D or 3D rotation matrix.');
end