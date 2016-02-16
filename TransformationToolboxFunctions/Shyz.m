function G = Shyz(S)
%Shyz homogeneous shear along the yz-direction
%   Shyz(S) creates a 3D shear transformation of magnitude S(1) along the
%   y-direction and S(2) along the z-direction.
%
%   See also Sx Sy Sz Shx Shy Shxy Shxz vee wedge
%
%   M. Kutzer 03Dec2014, USNA

G = [   1, 0, 0, 0;...
     S(1), 1, 0, 0;...
     S(2), 0, 1, 0;...
     0, 0,    0, 1];