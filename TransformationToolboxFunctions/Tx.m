function H = Tx(d,dim)
%Tx homogeneous translation along the z-axis
%   Tx(d) creates a 3D homogeneous transformation containing zero rotation 
%   and a translation of magnitude d along the x-axis.
%
%   Tx(d,3) creates a 3D homogeneous transformation containing zero  
%   rotation and a translation of magnitude d along the x-axis.
%
%   Tx(d,2) creates a 2D homogeneous transformation containing zero  
%   rotation and a translation of magnitude d along the x-axis.
%
%   See also Rx Ry Rz Ty Tz
%
%   M. Kutzer 20Oct2014, USNA

% Format enables use of symbolic variables

if nargin == 1
    dim = 3;
end

% 3D Homogeneous Transformation, SE(3)
if dim == 3
    H = [1, 0, 0, d;...
         0, 1, 0, 0;...
         0, 0, 1, 0;...
         0, 0, 0, 1];
    return
end

% 2D Homogeneous Transformation, SE(2)
if dim == 2
    H = [1, 0, d;...
         0, 1, 0;...
         0, 0, 1];
    return
end