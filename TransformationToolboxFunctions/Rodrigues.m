function R = Rodrigues(k,theta)
% RODRIGUES returns an element of the rotation group SO(3) given a
% direction and magnitude of rotation.
%   R = RODRIGUES(k,theta) defines the direction of rotation as unit vector
%   k and magnitude of rotation theta.
%
%   R = RODRIGUES(k) defines the direction and magnitude of rotation in a
%   single vector.
%
%   M. Kutzer, 28Jun2021, USNA

%% Check input(s)
narginchk(1,2)
if numel(k) ~= 3
    error('Axis of rotation must be defined as a 3-element vector.');
end

switch nargin
    case 1
        theta = norm(k);
        k = k./theta;
    case 2
        nk = norm(k);
        if nk ~= 1
            warning('Axis of rotation must be a unit vector, adjusting magnitude by %.11f to account for discrepency.',nk);
            theta = theta * nk;
            k = k./nk;
        end
end

%% Calculate rotation
K = wedge(k);
R = eye(3) + sin(theta)*K + (1-cos(theta))*(K^2);