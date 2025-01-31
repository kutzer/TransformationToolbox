function R = vectorToSO(v,dim)
% VECTORTOSO defines a rotation matrix given a vector
%   R = vectorToSO(v)
%
%   R = vectorToSO(v,dim)
%
%   Input(s)
%       v   - Nx1 array defining a unit vector direction in ND
%       dim - [OPTIONAL] positive integer specifying dimension of 
%             specified vector. Default dim = N.
%
%   Output(s)
%       R - NxN array defining rotation matrix with the column
%           corresponding to dim equal to v.
%
%           R(:,dim) = v
%
%   M. Kutzer, 31Jan2025, USNA

%% Check input(s)
narginchk(1,2);

N = numel(v);
if N < 2 || N > 3
    error('This function currently supports 2D and 3D vectors.');
end

if nargin < 2
    dim = N;
end

if dim ~= round(dim) || dim < 1 || dim > N
    error('Dimension must be specified as a positive integer between 1 and %d.',N);
end

%% Force v to unit vector
v_hat = reshape(v./norm(v),N,1);

%% Define rotation
switch N
    case 2
        w_hat = nCross(v_hat);

        switch dim
            case 1
                R = [+v_hat,-w_hat];
            case 2
                R = [+w_hat,+v_hat];
        end

    case 3
        % Define alternate direction
        w_hat = zeros(3,1);
        idx = find(abs(v_hat) == min(abs(v_hat)),1,'first');
        w_hat(idx) = 1;

        switch dim
            case 1
                x_hat = v_hat;
                y_est = w_hat;
                z_est = cross(x_hat,y_est);
                z_hat = z_est./norm(z_est);
                y_hat = cross(z_hat,x_hat);
            case 2
                y_hat = v_hat;
                z_est = w_hat;
                x_est = cross(y_hat,z_est);
                x_hat = x_est./norm(x_est);
                z_hat = cross(x_hat,y_hat);
            case 3
                z_hat = v_hat;
                x_est = w_hat;
                y_est = cross(z_hat,x_est);
                y_hat = y_est./norm(y_est);
                x_hat = cross(y_hat,z_hat);
        end
        R = [x_hat,y_hat,z_hat];
end
