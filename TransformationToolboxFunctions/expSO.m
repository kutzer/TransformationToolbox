function R = expSO(r)
% EXPSO calculates the matrix exponential of an element of the special
% orthogonal group.
%   R = EXPSO(r) calculates the matrix exponential using the general
%   formulation for SO(2), and Rodrigues's formula for SO(3). SO(N > 3)
%   uses the matrix exponential. 
%
%   Input(s)
%       r - NxN element of so(N)
%
%   Output(s)
%       R - NxN element of SO(N)
%
%   M. Kutzer, 26Jan2022, USNA

%% Convert to vector
k = vee(r);

%% Calculate axis/angle
Angle = norm(k);
if Angle == 0
    R = eye(size(r,1));
else
    Axis = transpose( k./Angle );
    R = AxisAngletoSO(Axis,Angle);
end