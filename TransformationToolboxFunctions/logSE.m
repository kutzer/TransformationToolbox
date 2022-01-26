function h = logSE(H)
% LOGSE calculates the matrix natural log of an element of the special
% Euclidean group.
%   h = LOGSE(H) calculates the matrix natural log using the general 
%   formulation for SE(2), and Rodrigues's formula for SE(3). SE(N > 3)
%   uses the matrix natural log.
%
%   Input(s)
%       H - (N+1)x(N+1) element of SE(N)
%
%   Output(s)
%       r - (N+1)x(N+1) element of se(N)
%
%   M. Kutzer, 26Jan2022, USNA

%% Default options
ZERO = 1e-8;

%% Check inputs
narginchk(1,1);
[bin,msg] = isSE(H,ZERO);
if ~bin
    warning('logSE:NotSE',...
        ['Input must be a valid n-dimensional rigid body transformation matrix.\n',...
        ' -> %s'],msg);
end

%% Calculate se(N)
N = size(H,1);
R = H(1:(N-1),1:(N-1));
d = H(1:(N-1),N);

h = zeros(N,N);
switch N
    case 3
        [Axis,Angle] = SOtoAxisAngle(R);
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
        [Axis,Angle] = SOtoAxisAngle(R);
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
        if isZero(R - eye(N-1))
            h(1:(N-1),N) = d;
        else
            h = logm(H);
        end
end