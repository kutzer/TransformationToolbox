function [H,H_i] = DHtableToFkin(DHtable)
% DHTABLETOFKIN calculates a transformation matrix (SE3) associated with 
% the forward kinematics defined in a DH table.
%   [H,H_i] = DHTABLETOFKIN(DHtable) This function creates a 4x4 array 
%   element of SE(3) representing the forward kinematics of a manipulator 
%   from a DH table following the Denavit–Hartenberg convention:
%
%   Input(s)
%       DHtable - Nx4 array containing the DH table parameters for the
%                 SIA20f for a given joint configuration q. 
%
%   Output(s)
%       H   - 4x4 element of SE(3) describing the pose of end-effector 
%             (or last frame) relative to base (or first frame).
%             H = H_{N}^{0}, see below
%       H_i - N-element cell array. Each element contains a 4x4 element of
%             SE(3) representing the relative transform between each
%             subsiquent frame.
%             H{i} = H_{i}^{i-1}, see below
%
%        DHtable = [theta_1,d_1,a_1,alpha_1; 
%                   theta_2,d_2,a_2,alpha_2;
%                   ...
%                   theta_N,d_N,a_N,alpha_N];
%
%        H_{i}^{i-1} = Rz(theta_i)*Tz(d_i)*Tx(a_i)*Rx(alpha_i);
%        H_{N}^{0} = H_{1}^{0}*H_{2}^{1}*...*H_{N}^{N-1};
%
%   See also DH plotDHtable
%   
%   M. Kutzer 14July2015, USNA

% Updates
%   07Jul2021 - Updated documentation and added error checking

%% Check inputs
narginchk(1,1);
N = size(DHtable,1);
M = size(DHtable,2);
if N < 1
    error('DH table must contain at least one row.');
end
if M ~= 4
    error('DH table must contain 4 columns.');
end

%% Calculate forward kinematics
H = eye(4);
N = size(DHtable,1);
for i = 1:N
    H_i{i} = DH(DHtable(i,1),DHtable(i,2),DHtable(i,3),DHtable(i,4));
	H = H*H_i{i};
end