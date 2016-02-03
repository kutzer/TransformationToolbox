function G = Sx(s,dim)
%Sx homogeneous scaling along the x-direction
%   Sx(s) creates a 3D scaling transformation of magnitude s along the
%   x-direction.
%
%   Sx(s,3) creates a 3D scaling transformation of magnitude s along the
%   x-direction.
%
%   Sx(s,2) creates a 2D scaling transformation of magnitude s along the
%   x-direction.
%
%   See also Sy Sz Shx Shy Shxy Shxz Shyz vee wedge
%
%   M. Kutzer 03Dec2014, USNA

if nargin == 1
    dim = 3;
end

if s == 0
    error('"s" should be non-zero.');
end

% 3D Scaling Transformation
if dim == 3
    G = [s, 0, 0, 0;...
         0, 1, 0, 0;...
         0, 0, 1, 0;...
         0, 0, 0, 1];
    return
end

% 2D Scaling Transformation
if dim == 2
    G = [s, 0, 0;...
         0, 1, 0;...
         0, 0, 1];
    return
end