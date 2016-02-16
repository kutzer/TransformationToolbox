function G = Shxy(S)
%Shxy homogeneous shear along the xy-direction
%   Shxy(S) creates a 3D shear transformation of magnitude S(1) along the
%   x-direction and S(2) along the y-direction.
%
%   See also Sx Sy Sz Shx Shy Shxz Shyz vee wedge
%
%   M. Kutzer 03Dec2014, USNA

G = [1, 0, S(1), 0;...
     0, 1, S(2), 0;...
     0, 0,    1, 0;...
     0, 0,    0, 1];