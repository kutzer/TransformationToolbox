function G = Shy(s)
%Shy homogeneous shear along the y-direction
%   SHy(s) creates a 2D shear transformation of magnitude s along the
%   y-direction.
%
%   See also Sx Sy Sz Shx Shxy Shxz Shyz vee wedge
%
%   M. Kutzer 03Dec2014, USNA

G = [1, 0, 0;...
     s, 1, 0;...
     0, 0, 1];
