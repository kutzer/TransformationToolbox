function r = logSO(R)
% LOGSO calculates the matrix natural log of an element of the special 
% orthogonal group.
%   r = LOGSO(R) calculates the matrix natural log using the general 
%   formulation for SO(2), and Rodrigues's formula for SO(3). SO(N > 3)
%   uses the matrix natural log.
%
%   Input(s)
%       R - NxN element of SO(N)
%
%   Output(s)
%       r - NxN element of so(N)
%
%   M. Kutzer 08Jan2016, USNA

% Updates
%   26Jan2022 - Replaced logm usage entirely
%   26Jan2022 - Updated documentation

%% Calculate logSO
[Axis,Angle] = SOtoAxisAngle(R);
v = Axis*Angle;
r = wedge(v);
