function G = Sz(s)
%Sz homogeneous scaling along the z-direction
%   Sz(s) creates a 3D scaling transformation of magnitude s along the
%   z-direction.
%
%   See also Sx Sy Shx Shy Shxy Shxz Shyz vee wedge
%
%   M. Kutzer 03Dec2014, USNA

% Updates
%   27Feb2016 - Updated help to correct typo

% 3D Scaling Transformation
G = [1, 0, 0, 0;...
     0, 1, 0, 0;...
     0, 0, s, 0;...
     0, 0, 0, 1];

