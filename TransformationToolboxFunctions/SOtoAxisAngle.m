function [Axis,Angle] = SOtoAxisAngle(R,varargin)
% SOtoAxisAngle converts a N-dimensional rotation matrix to to an 
% axis/angle parameterization.
%   [Axis,Angle] = SOtoAxisAngle(R)
%
%   [Axis,Angle] = SOtoAxisAngle(___,ZERO)
%
%   [Axis,Angle] = SOtoAxisAngle(___,fast)
%
%   Input(s)
%       R - NxN element of SO(N)
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not   
%              specified, a default value is used.
%       fast - [OPTIONAL] true/false logical value indicating whether to
%              skip checking SO(N). Choosing fast = true ignores specified
%              ZERO. 
%                fast = true    - Skip checking if R \in SO(N)
%                fast = [false] - Check if R \in SO(N)
%   Output(s)
%       Axis - scalar angle value in radians, Axis \in [0,\pi]
%
%   See also AxisAngletoSO logSO expSO logSE expSE
%
%   M. Kutzer 22Jan2016, USNA

% Updates
%   02Feb2016 - Updated to include n-dimensional axis/angle.
%   07Feb2018 - Updated to replace SO check with a warning.
%   26Jan2022 - Added increased "zero" value of 1e-8
%   26Jan2022 - Replaced "vrrotmat2vec" with "rotm2axang" for SO(3)
%   26Jan2022 - Added default axis for zero angle with SO(N) where N > 3
%   06Sep2022 - Updated to include ZERO and fast optional inputs
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast

% TODO - address negative eigenvalue issues of logm for larger than 3x3

%% Default options
ZERO = 1e-8;
fast = false;

%% Check inputs
narginchk(1,3);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms

%% Check R
if ~fast
    [bin,msg] = isSO(R,ZERO);
    if ~bin
        warning('SOtoAxisAngle:NotSO',...
            ['Input must be a valid n-dimensional rotation matrix.\n',...
            ' -> %s'],msg);
    end
end

%% Calculate axis angle
N = size(R,1);
switch N
    case 2
        Angle = atan2(R(2),R(1));
        Axis = sign(Angle);
        Angle = abs(Angle);
    case 3
        %r = vrrotmat2vec(R);
        r = rotm2axang(R);
        Axis = r(1:3);
        Angle = r(4);
    otherwise
        %error('SOtoAxisAngle:NotSO',...
        %    'Input must be a valid 2D or 3D rotation matrix.');
        
        % TODO - address negative eigenvalue issues of logm
        r = logm(R);
        v = vee(r,'fast');
        Angle = norm(v);
        if Angle ~= 0
            Axis = transpose( v./Angle );
        else
            Axis = zeros(1,numel(v));
            Axis(end) = 1;
        end
end