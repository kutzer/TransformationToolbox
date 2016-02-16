function H = Ry(theta)
%Ry homogeneous rotation about the y-axis (radians)
%   Ry(theta) creates a 3D homogeneous transformation containing a rotation
%   of magnitude theta about the y-axis and a translation of zero.
%   NOTE: "theta" must be defined in radians.
%
%   See also Rx Rz Tx Ty Tz vee wedge
%
%   M. Kutzer 20Oct2014, USNA

% Format enables use of symbolic variables
H(1,1) = cos(theta);
H(1,3) = sin(theta);
H(3,1) = -sin(theta);
H(3,3) = cos(theta);

H(2,2) = 1;
H(4,4) = 1;