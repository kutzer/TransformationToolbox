function DH = recoverDH(H)
% RECOVERDH calculate the parameters to populate a single row of a DH table
% given a 3D transformation. 
%
%   M. Kutzer 13Nov2014, USNA

ZERO = 1e-15;
%% Calculate z-rotation
theta = atan2(H(2,1),H(1,1));
alpha = atan2(H(3,2),H(3,3));

if ~strcmpi(class(theta),'sym')
    if abs(sin(theta)) > ZERO
        a = H(2,4)/sin(theta);
    else
        a = H(1,4)/cos(theta);
    end
else
    a = H(2,4)/sin(theta);
end

d = H(3,4);

DH = [theta,d,a,alpha];

%TODO - check the remaining matrix