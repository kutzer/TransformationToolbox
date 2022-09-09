function H = expSE(h,varargin)
% EXPSE calculates the matrix exponential of an element of the special
% Euclidean group.
%   H = EXPSE(h) calculates the matrix exponential using the general 
%   formulation for SE(2), and Rodrigues's formula for SE(3). SE(N > 3)
%   uses the matrix natural log.
%
%   H = EXPSE(___,ZERO)
%
%   H = EXPSE(___,FAST)
%
%   Input(s)
%       h      - (N+1)x(N+1) element of se(N)
%       ZERO   - [OPTIONAL] positive scalar value that is
%                sufficiently close to zero to be assumed zero
%                (e.g. ZERO = 1e-8). If a "ZERO" is not specified,
%                a default of ZERO = 1e-8 is used.
%       fast   - [OPTIONAL] true/false logical value indicating
%                whether to skip checking a specified property or
%                properties. If "fast" is not specified, a default of
%                fast = false is used.
%           fast = true    - Skip checking property or properties
%           fast = [false] - Check property or properties
%
%   Output(s)
%       r - (N+1)x(N+1) element of se(N)
%
%   M. Kutzer, 26Jan2022, USNA

% Update(s)
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast

%% Default options
ZERO = 1e-8;
fast = false;

%% Check inputs
narginchk(1,3);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

if ~fast
    % TODO - check h for valid se(N) properties
end

% TODO - Check bad inputs (cellOut)

%% Calculate se(N)
N = size(h,1);
r = h(1:(N-1),1:(N-1));
delta = h(1:(N-1),N);

% Convert to vector
k = vee(r,ZERO,fast);
% Calculate axis/angle
Angle = norm(k);
if Angle == 0
    R = eye(size(r,1));
else
    Axis = transpose( k./Angle );
    R = AxisAngletoSO(Axis,Angle);
end

H = eye(N);
H(1:(N-1),1:(N-1)) = R;
switch N
    case 3
        if Angle == 0
            H(1:(N-1),N) = delta;
        else
            K = wedge(Axis);
            J_L = eye((N-1)) +...
                ( (1-cos(Angle))/Angle )*K +...
                ( (Angle - sin(Angle))/Angle )*K^2;
            H(1:(N-1),N) = J_L * delta;
        end
    case 4
        if Angle == 0
            H(1:(N-1),N) = delta;
        else
            K = wedge(Axis);
            J_L = eye((N-1)) +...
                ( (1-cos(Angle))/Angle )*K +...
                ( (Angle - sin(Angle))/Angle )*K^2;
            H(1:(N-1),N) = J_L * delta;
        end
    otherwise
        if isZero(R - eye(N-1),ZERO)
            H(1:(N-1),N) = d;
        else
            % TODO - clean this up
            H = expm(h);
        end
end