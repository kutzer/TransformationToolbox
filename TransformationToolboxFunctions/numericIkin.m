function [q_Des,q_Init] = numericIkin(fkin,J_e,H_eDes2o,q_Init,s,delta_q_min)
% NUMERICIKIN solves the inverse kinematics problem numerically
%   q_Des = numericIkin(fkin,J_e,H_eDes2o,q_Init)
%
%   q_Des = numericIkin(fkin,J_e,H_eDes2o,q_Init,s)
%
%   q_Des = numericIkin(fkin,J_e,H_eDes2o,q_Init,[],delta_q_min)
%
%   q_Des = numericIkin(fkin,J_e,H_eDes2o,q_Init,s,delta_q_min)
%
%   Input(s)
%       fkin        - anonymous forward kinematics function
%       J_e         - anonymous body-fixed Jacobian function
%       H_eDes2o    - desired end-effector pose
%       q_Init      - Nx1 array describing initial joint configuration
%       s           - [OPTIONAL] scalar step size parameter (rad)
%       delta_q_min - [OPTIONAL] scalar minimum change in q (rad)
%
%   Output(s)
%       q_Des  - joint configuration associated with H_eDes2o
%       q_Init - joint configuration used to initialize inverse kinematics
%                solution. This value may differ from the supplied q_Now
%                value if a new 
%   M. Kutzer, 02Mar2022, USNA

debug = false;

%% Check input(s)
narginchk(4,6)

if nargin < 6
    delta_q_min = deg2rad(0.01);
end

if nargin < 5 || isempty(s)
    s = deg2rad(0.5);
end

% Impose Nx1 on q_Now
q_Init = reshape(q_Init,[],1);

% Test forward kinematics function
try
    H_e2o = fkin(q_Init);
catch ME
    error('The supplied forward kinematics function is not properly defined:\n\t"%s"',ME.message)
end

% Test inverse kinematics function
try
    J = J_e(q_Init);
catch ME
    error('The supplied forward kinematics function is not properly defined:\n\t"%s"',ME.message)
end

% Check Jacobian rows against number of task variables
n = size(H_e2o,1);
[M,N] = size(J);
if M ~= numel(seBasis(n-1))
    error('Jacobian is %dx%d, but the number of task parameters for SE(%d) is %d.',...
        M,N,n-1,numel(seBasis(n-1)));
end

% Check Jacobian columns against number of joint variables
if N ~= numel(q_Init)
    error('Jacobian is %dx%d, but the specified joint configuration is %d.',...
        M,N,numel(q_Init));
end

%% Move the arm to the desired configuration
% Define initial joint configuration
q_Now = q_Init;

delta_q_history = zeros(numel(q_Now),2);
while true
    % Update history with previous change in joint configuration
    delta_q_history(:,1) = delta_q_history(:,2);

    % Get current information from robot
    H_eNow2o = fkin(q_Now);    % current end-effector pose

    % Calculate change desired change in end-effector pose
    H_eDes2eNow = invSE(H_eNow2o)*H_eDes2o;

    % Calculate \Delta{\vec{v}}^e (relative to current end-effector frame)
    R_eDes2eNow = H_eDes2eNow(1:3,1:3);
    d_eDes2eNow = H_eDes2eNow(1:3,4);
    delta_k_eNow = vee( logSO(R_eDes2eNow),'fast' );
    delta_d_eNow = d_eDes2eNow;
    delta_v_eNow = [delta_k_eNow; delta_d_eNow];

    % Calculate delta_q using Jacobian
    delta_q = pinv( J_e(q_Now) ) * delta_v_eNow;

    % Scale delta_q
    inf_q = norm(delta_q,"inf");
    if inf_q > s
        delta_q = s*delta_q./inf_q;
    end

    % Update history with new change in joint configuration
    delta_q_history(:,2) = delta_q;

    % Check for oscillations
    % -> Maximum absolute value of change in delta_q
    inf_q_history = norm( sum(delta_q_history,2), "inf");
    % -> Difference in sign between past & present delta_q
    % NOTE: This value should be all zeros if the joint configuration is
    %       "wobbling"
    % (1) Round near-zero values to zero (for numeric stability)
    ZERO = (1e-1)*delta_q_min;
    tf_ZERO = abs(delta_q_history) < ZERO;
    delta_q_history(tf_ZERO) = 0;
    % (2) Check for alternating signs
    sgn_q_history = sum( sign(delta_q_history),2 );

    % Calculate updated q
    q_Now = q_Now + delta_q;

    % ---- Display debug information --------------------------------------
    if debug
        % Display delta_q
        fprintf('Change in joint config. ')
        fprintf('[')
        fprintf('%8.4f ',delta_q);
        fprintf('], |sum(hist)|_inf = %f\n',inf_q_history);
    end
    % ---------------------------------------------------------------------

    % Check for oscillations
    if inf_q_history < s && nnz(sgn_q_history) == 0
        % ---- Display debug information ----------------------------------
        if debug
            fprintf(2,'Oscillation detected, randomizing joint configuration\n')
            % Display current joint configuration
            fprintf(2,'\tCurrent joint config. ')
            fprintf('[')
            fprintf('%8.4f ',q_Now);
            fprintf(']\n');
        end
        % -----------------------------------------------------------------

        % Initialize new joint configuration
        % NOTE: There are more robust ways to choose a new initial joint
        %       configuration in this circumatance, but this should arrive
        %       at a solution eventually (if it exists)
        q_Now = 2*pi*rand(size(q_Now)) - pi;

        % Update initial joint configuration
        q_Init = q_Now;

        % ---- Display debug information ----------------------------------
        if debug
            % Display new joint configuration
            fprintf(2,'\t Random joint config. ')
            fprintf('[')
            fprintf('%8.4f ',q_Now);
            fprintf(']\n');
        end
        % -----------------------------------------------------------------
    end

    % Check if waypoint is achieved
    if norm( delta_q,"inf" ) < delta_q_min
        fprintf('Waypoint achieved\n');
        break
    end

end
q_Des = q_Now;