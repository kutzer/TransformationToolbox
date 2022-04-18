function H_e2o = fkin(q,L)
% UR3e forward kinematics function
%   H_e2o = fkin(q,L)
%
%   Input(s)
%       q - 6x1 array defining joint configuration (radians)
%       L - 6x1 array defining link lengths 
%
%   Output(s)
%       H_e2o - 4x4 array element of SE(3) defining the end-effector pose
%               of the UR3e relative to the base frame
%
%   M. Kutzer, 18Apr2022, USNA

% Published link lengths (reference)
% L = [151.85; 243.55; 213.20; 131.05; 85.35; 92.10];

H_1to0 = Rz(pi + q(1))*Rx(pi/2)*Ty(L(1));
H_2to1 = Rz(-q(2))*Tx(L(2));
H_3to2 = Rz(-q(3))*Tx(L(3));
H_4to3 = Rz(-q(4))*Rx(pi/2)*Ty(-L(4));
H_5to4 = Rz( q(5))*Rx(-pi/2)*Ty(-L(5));
H_6to5 = Rz(pi - q(6))*Rx(pi)*Tz(L(6));

H_e2o = H_1to0*H_2to1*H_3to2*H_4to3*H_5to4*H_6to5;