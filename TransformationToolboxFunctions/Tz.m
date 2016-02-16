function H = Tz(d)
%Tz homogeneous translation along the z-axis
%   Tz(d) creates a 3D homogeneous transformation containing zero rotation 
%   and a translation of magnitude d along the z-axis.
%
%   See also Rx Ry Rz Tx Ty
%
%   M. Kutzer 20Oct2014, USNA

% Format enables use of symbolic variables

H = [1, 0, 0, 0;...
     0, 1, 0, 0;...
     0, 0, 1, d;...
     0, 0, 0, 1];