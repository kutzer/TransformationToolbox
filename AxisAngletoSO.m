function R = AxisAngletoSO(Axis,Angle)
% AxisAngletoSO converts a 2D or 3D axis/angle combination to a rotation 
%   matrix (element of SO(2) or SO(3) to axis angle.
%
%   R = AxisAngletoSO(Axis,Angle)
%
%   M. Kutzer 22Jan2016, USNA

%% Check inputs
narginchk(2,2)
if numel(Axis) ~= 1 || numel(Axis) ~= 3
    error('AxisAngletoSO:BadAxis',...
        'An Axis must be specified as a 3-element array for 3D and a 1-element array for 2D.');
end
if numel(Angle) > 1
    error('AxisAngletoSO:BadAngle',...
        'The rotation angle must be specified as a scalar value.');
end

%% Calculate rotation
N = numel(Axis);
switch N
    case 1
        R = [cos(Axis*Angle), -sin(Axis*Angle); sin(Axis*Angle), cos(Axis*Angle)];
    case 3
        R = vrrotvec2mat([Axis,Angle]);
    otherwise
        error('AxisAngletoSO:BadAxis',...
            'An Axis must be specified as a 3-element array for 3D and a 1-element array for 2D.');
end