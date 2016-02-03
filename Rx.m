function H = Rx(theta)
%Rx homogeneous rotation about the x-axis (radians)
%   Rx(theta) creates a 3D homogeneous transformation containing a rotation
%   of magnitude theta about the x-axis and a translation of zero.
%   NOTE: "theta" must be defined in radians.
%
%   See also Ry Rz Tx Ty Tz vee wedge
%
%   M. Kutzer 20Oct2014, USNA

% Format enables use of symbolic variables
H(2,2) = cos(theta);
H(2,3) = -sin(theta);
H(3,2) = sin(theta);
H(3,3) = cos(theta);

H(1,1) = 1;
H(4,4) = 1;