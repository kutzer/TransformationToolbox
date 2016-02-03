function [H,H_i] = DHtableToFkin(DHtable)
% DHTABLETOFKIN calculates a transformation matrix (SE3) associated with 
% the forward kinematics defined in a DH table.
%   DHTABLETOFKIN(DHtable) This function creates a 4x4 array element of SE3
%   representing the forward kinematics of a manipulator from a DH table 
%   following the Denavit–Hartenberg convention:
%
%       DHtable = [theta_0,d_0,a_0,alpha_0; 
%                  theta_1,d_1,a_1,alpha_1;
%                  ...
%                  theta_N,d_N,a_N,alpha_N];
%       H_i = Rz(theta_i)*Tz(d_i)*Tx(a_i)*Rx(alpha_i);
%       H = H_0*H_1*...*H_N;
%
%   See also DH plotDHtable
%   
%   M. Kutzer 14July2015, USNA

H = eye(4);
N = size(DHtable,1);
for i = 1:N
    H_i{i} = DH(DHtable(i,1),DHtable(i,2),DHtable(i,3),DHtable(i,4));
	H = H*H_i{i};
end