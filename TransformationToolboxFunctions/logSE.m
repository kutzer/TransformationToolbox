function h = logSE(H,varargin)
% LOGSE calculates the matrix natural log of an element of the special
% Euclidean group.
%   h = LOGSE(H) calculates the matrix natural log using the general 
%   formulation for SE(2), and Rodrigues's formula for SE(3). SE(N > 3)
%   uses the matrix natural log.
%
%   h = LOGSE(___,ZERO)
%
%   h = LOGSE(___,fast)
%
%   Input(s)
%       H - (N+1)x(N+1) element of SE(N)
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not   
%              specified, a default value is used.
%       fast - [OPTIONAL] true/false logical value indicating whether to
%              skip checking SE(N). Choosing fast = true ignores specified
%              ZERO. 
%                fast = true    - Skip checking if H \in SE(N)
%                fast = [false] - Check if H \in SE(N)
%
%   Output(s)
%       r - (N+1)x(N+1) element of se(N)
%
%   See also expSE logSO expSO SOtoAxisAngle AxisAngletoSO 
%
%   M. Kutzer, 26Jan2022, USNA

% Updates:
%   06Sep2022 - Updated to include ZERO and fast optional inputs

%% Default options
ZERO = 1e-8;
fast = false;

%% Check inputs
narginchk(1,3);

if nargin > 1
    for i = 1:numel(varargin)
        switch lower( class(varargin{i} ))
            case 'logical'
                fast = varargin{i};
            otherwise
                if numel(varargin{i}) ~= 1
                    error('Numeric optional inputs must be scalar values.');
                end

                if varargin{i} == 0 || varargin{i} == 1
                    fast = logical(varargin{i});
                else
                    ZERO = varargin{i};
                end
        end
    end
end

% Check zero
if ZERO < 0
    error('ZERO value must be greater or equal to 0.')
end

%% Check H
if ~fast
    [bin,msg] = isSE(H,ZERO);
    if ~bin
        warning('logSE:NotSE',...
            ['Input must be a valid n-dimensional rigid body transformation matrix.\n',...
            ' -> %s'],msg);
    end
end

%% Calculate se(N)
N = size(H,1);
R = H(1:(N-1),1:(N-1));
d = H(1:(N-1),N);

h = zeros(N,N);
switch N
    case 3
        % Recover axis/angle
        [Axis,Angle] = SOtoAxisAngle(R);
        % Apply Rodrigues's Formula for SE(2)
        if Angle == 0
            h(1:(N-1),N) = d;
        else
            K = wedge(Axis);
            h(1:(N-1),1:(N-1)) = Angle*K;
            J_L = eye((N-1)) +...
                ( (1-cos(Angle))/Angle )*K +...
                ( (Angle - sin(Angle))/Angle )*K^2;
            invJ_L = (J_L)^-1;
            h(1:(N-1),N) = invJ_L * d;
        end
    case 4
        % Recover axis/angle
        [Axis,Angle] = SOtoAxisAngle(R);
        % Apply Rodrigues's Formula for SE(3)
        if Angle == 0
            h(1:(N-1),N) = d;
        else
            K = wedge(Axis);
            h(1:(N-1),1:(N-1)) = Angle*K;
            J_L = eye((N-1)) +...
                ( (1-cos(Angle))/Angle )*K +...
                ( (Angle - sin(Angle))/Angle )*K^2;
            invJ_L = (J_L)^-1;
            h(1:(N-1),N) = invJ_L * d;
        end
    otherwise
        % Apply generic SE(N) solution for N > 3 using logm
        if isZero(R - eye(N-1))
            h(1:(N-1),N) = d;
        else
            h = logm(H);
        end
end