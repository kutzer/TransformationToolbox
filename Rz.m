function H = Rz(theta,dim)
%Rz homogeneous rotation about the z-axis (radians)
%   Rz(theta) creates a 3D homogeneous transformation containing a rotation
%   of magnitude theta about the z-axis and a translation of zero. 
%   NOTE: "theta" must be defined in radians.
%
%   Rz(theta,3) creates a 3D homogeneous transformation containing a 
%   rotation of magnitude theta about the z-axis and a translation of zero. 
%   NOTE: "theta" must be defined in radians.
%
%   Rz(theta,2) creates a 2D homogeneous transformation containing a 
%   rotation of magnitude theta about the z-axis (assumed to be 
%   "out of the page") and a translation of zero. 
%   NOTE: "theta" must be defined in radians.
%
%   See also Rx Ry Tx Ty Tz vee wedge
%
%   M. Kutzer 20Oct2014, USNA

% Format enables use of symbolic variables
if nargin == 1
    dim = 3;
end

% 3D Homogeneous Transformation, SE(3)
if dim == 3
    H = [cos(theta), -sin(theta), 0, 0;...
         sin(theta),  cos(theta), 0, 0;...
                  0,           0, 1, 0;...
                  0,           0, 0, 1];
    return
end

% 2D Homogeneous Transformation, SE(2)
if dim == 2
    H = [cos(theta), -sin(theta), 0;...
         sin(theta),  cos(theta), 0;...
                  0,           0, 1];
    return
end