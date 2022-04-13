function [c,q,err] = approximateFixedKinematics(fcnH_e2o,H_e2o_ALL,options)
% APPROXIMATEFIXEDKINEMATICS approximates fixed values associated with an
% assumed forward kinematics function to best match a collected data set
% using fmincon.
%   [c,q,err] = approximateFixedKinematics(fcnH_e2o,H_e2o_ALL,options)
%
%   Input(s)
%       fcnH_e2o  - anonymous function defining forward kinematics.
%                   Function must use the following syntax:
%                       H_e2o = fcnH_e2o(q,c)
%       H_e2o_ALL - M-element cell array containing experimental data
%                   collected defining the pose of the "end-effector"
%                   (Frame e) relative to the "base" (Frame o).
%           H_e2o_ALL{i} - (N+1)x(N+1) array defining a valid element of
%                          SE(N)
%       options   - structured array containing options for approximating
%           *.qLims          - ax2 array containing lower/upper bounds of
%                              joint values
%           *.q0             - ax1 array containing initial joint values
%           *.cLims          - bx2 array containing lower/upper bounds of
%                              kinematic constants
%           *.c0             - bx1 array containing initial values of
%                              kinematic constants
%           *.RotationWeight - positive scalar term defining weight for
%                              rotation/orientation alignment
%
%   Output(s)
%       c   - nx1 array containing fixed values of fcnH_e2o recovered from
%             the data set.
%       q   - mxM array containing the joint configurations associated with
%             fcnH_e2o recovered from the data set.
%       err - structured array containing error statistics of fixed
%             kinematics.
%           *.Mean      - mean error 
%           *.Max       - max error
%           *.Min       - min error
%           *.STD       - error standard deviation
%           *.Veriance  - error variance
%
%       NOTE: q(:,i) corresponds to H_e2o_ALL{i} such that
%               H_e2o_ALL{i} \approx fcnH_e2o(q(i,:),c)
%
%       err - scalar value describing mean error
%
%   M. Kutzer, 12Apr2022, USNA

%% Check input(s)
narginchk(3,3);
fields = {'qLims','q0','cLims','c0','RotationWeight'};
bin = isfield(options,fields);
if nnz(bin) ~= numel(fields)
    str = sprintf('"options" is missing the following fields:\n');
    for i = reshape(find(bin),1,[])
        str = sprintf('%s\t%14s',str,fields{i});
    end
end

% Check for correct sizes
a = numel(options.q0);
if size(options.qLims,1) ~= a || size(options.qLims,2) ~= 2
    error('"options.qLims" must be %dx2 given a %dx1 value for "options.q0"',...
        a,a);
end

b = numel(options.c0);
if size(options.cLims,1) ~= b || size(options.cLims,2) ~= 2
    error('"options.cLims" must be %dx2 given a %dx1 value for "options.c0"',...
        a,a);
end

% Check forward kinematics function
q = rand(a,1);
c = rand(b,1);
try
    H_e2o = fcnH_e2o(q,c);
catch
    error('"fcnH_e2o" must be defined using the following syntax:\n\tH_e2o = fcnH_e2o(q,c);');
end
% TODO - check H_e2o against values of H_e2o_ALL

% TODO - check remaining input(s)

%% Define symbolic forward kinematics
% Define symbolic joint variables
q = sym('q',[a,1]);
% Define symbolic constants
c = sym('c',[b,1]);
% Define symbolic forward kinematics
H_e2o_sym = fcnH_e2o(q,c);

%% Calculate Jacobians for gradient definitions
J_qo = calculateJacobian(q,H_e2o_sym,'Constants',c,...
    'Order','RotationTranslation','Reference','Base');
J_co = calculateJacobian(c,H_e2o_sym,'Constants',q,...
    'Order','RotationTranslation','Reference','Base');

%% Setup Optimization Parameters
% Optimization options
fOptions = optimoptions(@fminunc,'Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true,'HessianFcn','objective');

% Optimization parameters
N = numel(H_e2o_ALL);
q0   = repmat(options.q0,1,N);
q_lb = repmat(options.qLims(:,1),1,N);
q_ub = repmat(options.qLims(:,2),1,N);
c0   = options.c0;
c_lb = options.cLims(:,1);
c_ub = options.cLims(:,2);

% Reshape optimization parameters
q0   = reshape(q0  ,1,[]);
q_lb = reshape(q_lb,1,[]);
q_ub = reshape(q_ub,1,[]);
c0   = reshape(c0  ,1,[]);
c_lb = reshape(c_lb,1,[]);
c_ub = reshape(c_ub,1,[]);

% Optimization functions
%fcost_q = @(q_ALL)cost_q(q_ALL,a,c,fcnH_e2o,J_qo,H_e2o_ALL,options.RotationWeight);
%fcost_c = @(c)cost_c(q_ALL,a,c,fcnH_e2o,J_co,H_e2o_ALL,options.RotationWeight);

%% Optimize
% fmincon - constrained minimization
% fminunc - unconstrained minimization
q = q0;
c = c0;

% Loop through q-minimization/c-minimization until... TBD stop condition
iter = 0;
while true
    fcost_q = @(q_ALL)cost_q(q_ALL,a,c,fcnH_e2o,J_qo,H_e2o_ALL,options.RotationWeight);
    [q,fval,exitflag,output] = fmincon(fcost_q,q,[],[],[],[],q_lb,q_ub,[],options);
    
    fcost_c = @(c)cost_c(q,a,c,fcnH_e2o,J_co,H_e2o_ALL,options.RotationWeight);
    [c,fval,exitflag,output] = fmincon(fcost_c,c,[],[],[],[],c_lb,c_ub,[],options);
    iter = iter+1;

    if iter == 3
        break
    end
end

%% Package outputs
[err,q,c] = errStats(q,a,c,fcnH_e2o,H_e2o_ALL,options.RotationWeight);

end

%% Define cost functions
function [cost,costEq,grad,gradEq] = cost_q(q_ALL,a,c,fcnH_e2o,J_qo,H_e2o_ALL,rw)
% Reshape input(s)
q_ALL = reshape(q_ALL,a,[]);
c = reshape(c,[],1);
% Calculate forward kinematics & estimate error
n = size(q_ALL,2);
cost = 0;
costEq = [];
if nargout > 2
    grad = zeros( size(q_ALL) );
    gradEq = [];
end
idx = 0;
for i = 1:n
    % --- Define Cost ---
    % Recover "true" (experimental) value
    H_eTru2o = H_e2o_ALL{i};
    H_o2eTru = invSE(H_eTru2o);
    % Recover "estimated" value
    H_eEst2o = fcnH_e2o(q_ALL(:,i),c);
    % Estimate error
    H_eEst2eTru = H_o2eTru*H_eEst2o;
    % Recover rotation & translation errors
    R_eEst2eTru = H_eEst2eTru(1:(end-1),1:(end-1));
    d_eEst2eTru = H_eEst2eTru(1:(end-1),end);
    rotErr = rw*R_eEst2eTru;
    trnErr = d_eEst2eTru;
    % Combine errors
    err = rotErr + repmat(trnErr,1,size(rotErr,2));
    err = diag( err*err.' );
    cost = cost + sum(err);

    % --- Define Gradient ---
    if nargout > 2
        J = J_qo(q_ALL(:,i),c);
        for k = 1:size(J,2)
            v = J(:,k);
            % TODO - expand to include SE(N), current is SE(3)
            % Derivative term (referenced to base)
            dR_eEst2o = expSO( wedge(v(1:3,1)));
            dd_eEst2o = v(4:6,1);
            dH_eEst2o = zeros(4);
            dH_eEst2o(1:3,1:3) = dR_eEst2o;
            dH_eEst2o(1:3,4) = dd_eEst2o;
            % Derivative term (referenced to eTru)
            dH_eEst2eTru = H_o2eTru*dH_eEst2o;
            % Recover rotation & translation gradient errors
            dR_eEst2eTru = dH_eEst2eTru(1:(end-1),1:(end-1));
            dd_eEst2eTru = dH_eEst2eTru(1:(end-1),end);
            dRotErr = rw*dR_eEst2eTru;
            dTrnErr = dd_eEst2eTru;
            % Combine gradient errors
            dErr = dRotErr + repmat(dTrnErr,1,size(dRotErr,2));
            idx = idx+1;
            grad(idx) = diag( dErr*err.' + err*dErr.' )/n;
        end
    end
end
cost = cost/n;
end

function [cost,costEq,grad,gradEq] = cost_c(q_ALL,a,c,fcnH_e2o,J_co,H_e2o_ALL)
% Reshape input(s)
q_ALL = reshape(q_ALL,a,[]);
c = reshape(c,[],1);
% Calculate forward kinematics & estimate error
n = size(q_ALL,2);
cost = 0;
costEq = [];
if nargout > 2
    grad = zeros( size(c) );
    gradEq = [];
end
%idx = 0;
for i = 1:n
    % --- Define Cost ---
    % Recover "true" (experimental) value
    H_eTru2o = H_e2o_ALL{i};
    H_o2eTru = invSE(H_eTru2o);
    % Recover "estimated" value
    H_eEst2o = fcnH_e2o(q_ALL(:,i),c);
    % Estimate error
    H_eEst2eTru = H_o2eTru*H_eEst2o;
    % Recover rotation & translation errors
    R_eEst2eTru = H_eEst2eTru(1:(end-1),1:(end-1));
    d_eEst2eTru = H_eEst2eTru(1:(end-1),end);
    rotErr = rw*R_eEst2eTru;
    trnErr = d_eEst2eTru;
    % Combine errors
    err = rotErr + repmat(trnErr,1,size(rotErr,2));
    err = diag( err*err.' );
    cost = cost + sum(err);

    % --- Define Gradient ---
    if nargout > 2
        J = J_co(c,q_ALL(:,i));
        for k = 1:size(J,2)
            v = J(:,k);
            % TODO - expand to include SE(N), current is SE(3)
            % Derivative term (referenced to base)
            dR_eEst2o = expSO( wedge(v(1:3,1)));
            dd_eEst2o = v(4:6,1);
            dH_eEst2o = zeros(4);
            dH_eEst2o(1:3,1:3) = dR_eEst2o;
            dH_eEst2o(1:3,4) = dd_eEst2o;
            % Derivative term (referenced to eTru)
            dH_eEst2eTru = H_o2eTru*dH_eEst2o;
            % Recover rotation & translation gradient errors
            dR_eEst2eTru = dH_eEst2eTru(1:(end-1),1:(end-1));
            dd_eEst2eTru = dH_eEst2eTru(1:(end-1),end);
            dRotErr = rw*dR_eEst2eTru;
            dTrnErr = dd_eEst2eTru;
            % Combine gradient errors
            dErr = dRotErr + repmat(dTrnErr,1,size(dRotErr,2));
            grad(k) = grad(k) + diag( dErr*err.' + err*dErr.' )/n;
        end
    end
end
cost = cost/n;
end

function [err,q_ALL,c] = errStats(q_ALL,a,c,fcnH_e2o,H_e2o_ALL,rw)
% Reshape input(s)
q_ALL = reshape(q_ALL,a,[]);
c = reshape(c,[],1);
% Calculate forward kinematics & estimate error
n = size(q_ALL,2);
err_ALL = zeros(1,n);
for i = 1:n
    % --- Define Cost ---
    % Recover "true" (experimental) value
    H_eTru2o = H_e2o_ALL{i};
    H_o2eTru = invSE(H_eTru2o);
    % Recover "estimated" value
    H_eEst2o = fcnH_e2o(q_ALL(:,i),c);
    % Estimate error
    H_eEst2eTru = H_o2eTru*H_eEst2o;
    % Recover rotation & translation errors
    R_eEst2eTru = H_eEst2eTru(1:(end-1),1:(end-1));
    d_eEst2eTru = H_eEst2eTru(1:(end-1),end);
    rotErr = rw*R_eEst2eTru;
    trnErr = d_eEst2eTru;
    % Combine errors
    err = rotErr + repmat(trnErr,1,size(rotErr,2));
    err = diag( err*err.' );
    err_ALL = sum(err);
end

err.Mean = mean(err_ALL);
err.Max = max(err_ALL);
err.Min = min(err_ALL);
err.STD  = std(err_ALL);
err.Variance = var(err_ALL);

end
