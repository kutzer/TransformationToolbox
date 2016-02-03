function G = Shxz(S)
%Shxz homogeneous shear along the xz-direction
%   Shxz(S) creates a 3D shear transformation of magnitude S(1) along the
%   x-direction and S(2) along the z-direction.
%
%   See also Sx Sy Sz Shx Shy Shxy Shyz vee wedge
%
%   M. Kutzer 03Dec2014, USNA

G = [1, S(1), 0, 0;...
     0,    1, 0, 0;...
     0, S(2), 1, 0;...
     0, 0,    0, 1];