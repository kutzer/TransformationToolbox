function DH = recoverDH(H,ZERO)
% RECOVERDH calculate the parameters to populate a single row of a DH table
% given a 3D transformation. 
%   DH = RECOVERDH(H)
%   DH = RECOVERDH(H,ZERO)
%
%   Input(s)
%       H  - 4x4 array element of SE(3) 
%       ZERO - Optional input allowing user to specify value a value that
%              is close enough to 0 to be considered 0. Default value is
%              1e-5. 
%   Output(s)
%       DH - 1x4 array containing a single row of a DH table such that 
%            DH = [theta, d, a, alpha]
%
%   Note: H = Rz(theta)*Tz(d)*Tx(a)*Rx(alpha). If the provided
%         transformation cannot be parameterized using DH, and empty set 
%         is returned
%
%   M. Kutzer 13Nov2014, USNA

% Update(s)
%   29Jun2021 - Updated documentation
%   29Jun2021 - Returns an empty set if the result is not transform that
%               can be represented using DH

%% Check inputs
narginchk(1,2);
% TODO - check if H \in SE (note isSE does not support symbolic)

if nargin > 1
    if numel(ZERO) ~= 1 || ZERO(1) < 0
        error('ZERO parameter must be a scalar value greater than 0.');
    end
end

%% Set default(s)
if nargin < 2
    ZERO = 1e-5;
end

%% Check matrix for DH form:
% [cos(theta),-cos(alpha)*sin(theta), sin(alpha)*sin(theta),a*cos(theta)]
% [sin(theta), cos(alpha)*cos(theta),-sin(alpha)*cos(theta),a*sin(theta)]
% [         0,            sin(alpha),            cos(alpha),           d]
% [         0,                     0,                     0,           1]
warn_str = 'Transformation does not appear to be associated with a single row of a DH table.';

chk = zeroFPError( H(2,1).^2 + H(1,1).^2, ZERO );
if chk ~= 1
    warning(warn_str);
    DH = [];
    return
end
    
chk = zeroFPError( H(3,2).^2 + H(3,3).^2, ZERO );
if chk ~= 1
    warning(warn_str);
    DH = [];
    return
end

aa1 = zeroFPError( H(1,4)./H(1,1), ZERO );
aa2 = zeroFPError( H(2,4)./H(2,1), ZERO );
if abs( zeroFPError(aa1 - aa2, ZERO) ) > 0 
    warning(warn_str);
    DH = [];
    return
end

%% Calculate z-rotation
theta = atan2(H(2,1),H(1,1));
alpha = atan2(H(3,2),H(3,3));

if ~strcmpi(class(theta),'sym')
    %if abs(sin(theta)) > abs(cos(theta))
    %    a = H(2,4)/sin(theta);
    %else
    %    a = H(1,4)/cos(theta);
    %end
    a = sqrt( H(2,4).^2 + H(1,4).^2 );
else
    %a = H(2,4)/sin(theta);
    a = sqrt( H(2,4).^2 + H(1,4).^2 );
end

d = H(3,4);

DH = [theta,d,a,alpha];
DH = zeroFPError(DH,ZERO);
%% Check result
% - atan2( sin(theta), cos(theta) ) will not simplify!!!
% H_chk = Rz(theta)*Tz(d)*Tx(a)*Rx(alpha);
% chk = simplify( real( H - H_chk ))
% chk = zeroFPError(H - H_chk)