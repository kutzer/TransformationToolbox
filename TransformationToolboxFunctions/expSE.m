function H = expSE(h)
% EXPSE calculates the matrix exponential of an element of the special
% Euclidean group.
%   H = EXPSE(h) calculates the matrix exponential using the general 
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
% TODO - Check inputs

%% Calculate se(N)
N = size(h,1);
r = h(1:(N-1),1:(N-1));
delta = h(1:(N-1),N);

% Convert to vector
k = vee(r);
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
        if isZero(R - eye(N-1))
            H(1:(N-1),N) = d;
        else
            % TODO - clean this up
            H = expm(h);
        end
end