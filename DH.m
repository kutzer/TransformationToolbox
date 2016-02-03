function H = DH(theta,d,a,alpha)
% DH calculates a transformation matrix (SE3) from four DH parameters.
%   DH(theta,d,a,alpha) This function creates a 4x4 array element of SE3
%   from four parameters following Denavit–Hartenberg convention:
%       H = Rz(theta)*Tz(d)*Tx(a)*Rx(alpha)
%
%   See also DHtableToFkin
%   
%   M. Kutzer 12Nov2014, USNA

H = Rz(theta)*Tz(d)*Tx(a)*Rx(alpha);