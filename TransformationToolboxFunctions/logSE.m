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
%              specified, a default value of ZERO = 1e-8 is used.
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
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast
%   28Nov2023 - Updated to include additional notes on left Jacobian

%% Default options
ZERO = 1e-8;
fast = false;

%% Check inputs
narginchk(1,3);

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms

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
    case {3,4}
        % Recover axis/angle
        [Axis,Angle] = SOtoAxisAngle(R,ZERO,fast);
        % Apply Rodrigues's Formula for SE(2) & SE(3)
        if Angle == 0
            h(1:(N-1),N) = d;
        else
            K = wedge(Axis);
            h(1:(N-1),1:(N-1)) = Angle*K;
            % Original method
            J_L = eye((N-1)) +...
                ( (1-cos(Angle))/Angle )*K +...
                ( (Angle - sin(Angle))/Angle )*K^2;
            % Alternative method
            %   Note: 
            %       Using this method does not provide the same result
            %       logm! 
            %
            %   References:
            %       [1] https://www.researchgate.net/publication/261711710_Computing_the_Logarithm_of_Homogenous_Matrices_in_SE3
            %       [2] https://natanaso.github.io/ece276a2019/ref/ECE276A_12_SE3.pdf
            %{
            J_L = eye((N-1)) +...
                ( (1-cos(Angle))/(Angle^2) )*K +...
                ( (Angle - sin(Angle))/(Angle^3) )*K^2;
            %}

            invJ_L = (J_L)^-1;
            h(1:(N-1),N) = invJ_L * d;
        end
    otherwise
        % Apply generic SE(N) solution for N > 3 using logm
        if isZero(R - eye(N-1),ZERO)
            h(1:(N-1),N) = d;
        else
            h = logm(H);
        end
end