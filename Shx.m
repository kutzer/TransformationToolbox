function G = Shx(s)
%Shx homogeneous shear along the x-direction
%   Shx(s) creates a 2D shear transformation of magnitude s along the
%   x-direction.
%
%   Shx(s,2) creates a 2D shear transformation of magnitude s along the
%   x-direction.
%
%   See also Sx Sy Sz Shy Shxy Shxz Shyz vee wedge
%
%   M. Kutzer 03Dec2014, USNA

G = [1, s, 0;...
     0, 1, 0;...
     0, 0, 1];
