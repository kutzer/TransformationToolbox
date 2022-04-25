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
%           *.q0             - axM array containing initial joint values
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

% Updates
%   25Apr2022 - q0 is now axM for M-element H_e2o_ALL

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
M = numel(H_e2o_ALL);
a = size(options.q0,1);
if M ~= size(options.q0,2)
    error('"options.q0" must contain %d columns given %d elements of "H_e2o_ALL',...
        M,M);
end

if size(options.qLims,1) ~= a || size(options.qLims,2) ~= 2
    error('"options.qLims" must be %dx2 given a %dx%d value for "options.q0"',...
        a,a,M);
end

b = numel(options.c0);
if size(options.cLims,1) ~= b || size(options.cLims,2) ~= 2
    error('"options.cLims" must be %dx2 given a %dx1 value for "options.c0"',...
        b,b);
end

% Check forward kinematics function
q = rand(a,1);
c = rand(b,1);
try
    H_e2o = fcnH_e2o(q,c);
catch
    error(...
        sprintf('"fcnH_e2o" must be defined using the following syntax:\n\tH_e2o = fcnH_e2o(q,c);')...
        );
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
J_qco = calculateJacobian([q; c],H_e2o_sym,...
    'Order','RotationTranslation','Reference','Base');

%% Setup Optimization Parameters
% Optimization options
%fOptions = optimoptions(@fminunc,'Algorithm','trust-region',...
%    'SpecifyObjectiveGradient',true,'HessianFcn','objective');
% fOptions_q = optimoptions(@fmincon,'Algorithm','trust-region-reflective',...    
%     'SpecifyObjectiveGradient',true,'CheckGradients',false,...
%     'Display','off',...
%     'PlotFcn',@(x,optimValues,state)statusPlot(x,optimValues,state,'q'));
fOptions_q = optimoptions(@fmincon,'Algorithm','interior-point',...    
    'Display','off',...
    'PlotFcn',@(x,optimValues,state)statusPlot(x,optimValues,state,'q'));
fOptions_c = optimoptions(@fmincon,'Algorithm','trust-region-reflective',...
    'SpecifyObjectiveGradient',true,'CheckGradients',false,...
    'Display','off',...
    'PlotFcn',@(x,optimValues,state)statusPlot(x,optimValues,state,'c'));
fOptions_qc = optimoptions(@fmincon,'Algorithm','interior-point',...    
    'Display','off',...
    'PlotFcn',@(x,optimValues,state)statusPlot(x,optimValues,state,'qc'));

% Optimization parameters
%M = numel(H_e2o_ALL);
%q0   = repmat(options.q0,1,M); % TODO: consider implementing if q0 has 1 column
q0 = options.q0;
q_lb = repmat(options.qLims(:,1),1,M);
q_ub = repmat(options.qLims(:,2),1,M);
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

% Define combine parameters
qc0 = [q0, c0];
qc_lb = [q_lb, c_lb];
qc_ub = [q_ub, c_ub];

% Optimization functions
%fcost_q = @(q_ALL)cost_q(q_ALL,a,c,fcnH_e2o,J_qo,H_e2o_ALL,options.RotationWeight);
%fcost_c = @(c)cost_c(q_ALL,a,c,fcnH_e2o,J_co,H_e2o_ALL,options.RotationWeight);

%% Optimize
% fmincon - constrained minimization
% fminunc - unconstrained minimization
q = q0;
c = c0;

% Loop through q-minimization/c-minimization until... TBD stop condition
fig = figure('Name','Optimization');
axs = axes('Parent',fig);
hold(axs,'on');
plt = plot(axs,nan,nan,'k');
plt_q = plot(axs,nan,nan,'.g');
plt_c = plot(axs,nan,nan,'.b');
plt_qc = plot(axs,nan,nan,'.r');
x = [];
y = [];
y_q = [];
y_c = [];
y_qc = [];

iter = 0;
while true
    % --- q ---
    fcost_q = @(q)cost_q(q,a,c,fcnH_e2o,J_qo,H_e2o_ALL,options.RotationWeight);
    [q,fval,exitflag,output] = fmincon(fcost_q,q,[],[],[],[],q_lb,q_ub,[],fOptions_q);
    iter = iter + (1/3);
    x(end+1) = iter;
    y(end+1) = fval;
    y_q(end+1) = fval;
    y_c(end+1) = nan;
    y_qc(end+1) = nan;

    set(plt,'XData',x,'YData',y);
    set(plt_q,'XData',x,'YData',y_q);
    set(plt_c,'XData',x,'YData',y_c);
    set(plt_qc,'XData',x,'YData',y_qc);

    % --- c ---
    fcost_c = @(c)cost_c(q,a,c,fcnH_e2o,J_co,H_e2o_ALL,options.RotationWeight);
    [c,fval,exitflag,output] = fmincon(fcost_c,c,[],[],[],[],c_lb,c_ub,[],fOptions_c);
    iter = iter + (1/3);
    x(end+1) = iter;
    y(end+1) = fval;
    y_q(end+1) = nan;
    y_c(end+1) = fval;
    y_qc(end+1) = nan;

    set(plt,'XData',x,'YData',y);
    set(plt_q,'XData',x,'YData',y_q);
    set(plt_c,'XData',x,'YData',y_c);
    set(plt_qc,'XData',x,'YData',y_qc);

    % --- qc ---
    qc = [q,c];
    fcost_qc = @(qc)cost_qc(qc,a,b,fcnH_e2o,J_qco,H_e2o_ALL,options.RotationWeight);
    [qc,fval,exitflag,output] = fmincon(fcost_qc,qc,[],[],[],[],qc_lb,qc_ub,[],fOptions_qc);
    iter = iter + (1/3);
    x(end+1) = iter;
    y(end+1) = fval;
    y_q(end+1) = nan;
    y_c(end+1) = nan;
    y_qc(end+1) = fval;

    set(plt,'XData',x,'YData',y);
    set(plt_q,'XData',x,'YData',y_q);
    set(plt_c,'XData',x,'YData',y_c);
    set(plt_qc,'XData',x,'YData',y_qc);
    
    c = qc( (end-b+1:end) );
    q = qc( 1:(end-b) );


    if iter == 1000
        break
    end
end
obj.Figure = findobj('Type','figure','Name','Optimization Plot Function');
delete(obj.Figure);

%% Package outputs
[err,q,c] = errStats(q,a,c,fcnH_e2o,H_e2o_ALL,options.RotationWeight);

end

%% Define cost functions
%function [cost,costEq,grad,gradEq] = cost_q(q_ALL,a,c,fcnH_e2o,J_qo,H_e2o_ALL,rw)
function [cost,grad] = cost_q(q_ALL,a,c,fcnH_e2o,J_qo,H_e2o_ALL,rw)
% Reshape input(s)
q_ALL = reshape(q_ALL,a,[]);
c = reshape(c,[],1);
% Calculate forward kinematics & estimate error
n = size(q_ALL,2);
cost = 0;
%costEq = [];
%if nargout > 2
if nargout > 1
    grad = zeros( size(q_ALL) );
    %gradEq = [];
end

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
    %if nargout > 2
    if nargout > 1
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

            grad(k,i) =...
                grad(k,i) +...
                sum( diag( dErr*diag(err).' + diag(err)*dErr.' ) )/n;
        end
    end
end
cost = cost/n;

% Reshape gradient for fmincon
%if nargout > 2
if nargout > 1
    grad = reshape(grad,1,[]);
end

end

%function [cost,costEq,grad,gradEq] = cost_c(q_ALL,a,c,fcnH_e2o,J_co,H_e2o_ALL,rw)
function [cost,grad] = cost_c(q_ALL,a,c,fcnH_e2o,J_co,H_e2o_ALL,rw)
% Reshape input(s)
q_ALL = reshape(q_ALL,a,[]);
c = reshape(c,[],1);
% Calculate forward kinematics & estimate error
n = size(q_ALL,2);
cost = 0;
%costEq = [];
%if nargout > 2
if nargout > 1
    grad = zeros( size(c) );
    %gradEq = [];
end

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
            grad(k) =...
                grad(k) +...
                sum( diag( dErr*diag(err).' + diag(err)*dErr.' ) )/n;
        end
    end
end
cost = cost/n;
end


function [cost,grad] = cost_qc(qc_ALL,a,b,fcnH_e2o,J_qco,H_e2o_ALL,rw)
% Parse input
c = qc_ALL( (end-b+1:end) );
q_ALL = qc_ALL( 1:(end-b) );
% Reshape inputs
q_ALL = reshape(q_ALL,a,[]);
c = reshape(c,[],1);
% Calculate forward kinematics & estimate error
n = size(q_ALL,2);
cost = 0;
%costEq = [];
%if nargout > 2
if nargout > 1
    grad = zeros( size(qc_ALL) );
    %gradEq = [];
end

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
    %if nargout > 2
    if nargout > 1
        J = J_qco([q_ALL(:,i); c]);
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

            grad(k,i) =...
                grad(k,i) +...
                sum( diag( dErr*diag(err).' + diag(err)*dErr.' ) )/n;
        end
    end
end
cost = cost/n;

% Reshape gradient for fmincon
%if nargout > 2
if nargout > 1
    grad = reshape(grad,1,[]);
end

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

clearvars err
err.Mean = mean(err_ALL);
err.Max = max(err_ALL);
err.Min = min(err_ALL);
err.STD  = std(err_ALL);
err.Variance = var(err_ALL);

end



function stop = statusPlot(~,optimValues,state,method)%varargin)

persistent obj

% optimValues
% state
% for i = 1:numel(varargin)
%     varargin{i}
% end

stop = false;

% Populate obj
if isempty(obj)
    % Create figure, axes, etc
    %obj.Figure = figure('Name','Optimization Cost','NumberTitle','off',...
    %    'ToolBar','none','MenuBar','None');
    %obj.Axes = axes('Parent',obj.Figure);
    
    obj.Figure = findobj('Type','figure','Name','Optimization Plot Function');
    obj.Axes = findobj('Type','axes','Parent',obj.Figure);
    hold(obj.Axes,'on');
    obj.Line_q = plot(obj.Axes,nan,nan,'g','tag','disable');
    obj.Line_c = plot(obj.Axes,nan,nan,'b','tag','disable');
end

% Populate objects
if ~ishandle(obj.Figure)
    obj.Figure = findobj('Type','figure','Name','Optimization Plot Function');
    obj.Axes = findobj('Type','axes','Parent',obj.Figure);
    hold(obj.Axes,'on');
    obj.Line_q = plot(obj.Axes,nan,nan,'g','tag','disable');
    obj.Line_c = plot(obj.Axes,nan,nan,'b','tag','disable');
end

if ~ishandle(obj.Axes)
    obj.Axes = findobj('Type','axes','Parent',obj.Figure);
    hold(obj.Axes,'on');
    obj.Line_q = plot(obj.Axes,nan,nan,'g','tag','disable');
    obj.Line_c = plot(obj.Axes,nan,nan,'b','tag','disable');
end

if ~ishandle(obj.Line_q)
    obj.Line_q = plot(obj.Axes,nan,nan,'g','tag','disable');
end

if ~ishandle(obj.Line_c)
    obj.Line_c = plot(obj.Axes,nan,nan,'b','tag','disable');
end

switch state
    case 'iter'
        switch lower(method)
            case 'q'
                x = get(obj.Line_q,'XData');
                y = get(obj.Line_q,'YData');
                
                switch lower( get(obj.Line_q,'Tag') )
                    case 'disable'
                        x(end+1) = nan;
                        y(end+1) = nan;
                    case 'enable'
                        % plot is already enabled
                    otherwise
                        error('Unexpected tag: %s',get(obj.Line_q,'Tag'));
                end
                set(obj.Line_q,'Tag','enable');
                set(obj.Line_c,'Tag','disable');

                i0 = find(~isnan(x),1,'last');
                if ~isempty(i0)
                    x0 = x(i0);
                else
                    x0 = 0;
                end

                x(end+1) = x0 + optimValues.iteration;
                y(end+1) = optimValues.fval;

                set(obj.Line_q,'XData',x,'YData',y);
            %case 'c'
            otherwise
                x = get(obj.Line_c,'XData');
                y = get(obj.Line_c,'YData');

                switch lower( get(obj.Line_c,'Tag') )
                    case 'disable'
                        x(end+1) = nan;
                        y(end+1) = nan;

                    case 'enable'
                        % plot is already enabled
                    otherwise
                        error('Unexpected tag: %s',get(obj.Line_q,'Tag'));
                end

                set(obj.Line_q,'Tag','disable');
                set(obj.Line_c,'Tag','enable');

                i0 = find(~isnan(x),1,'last');
                if ~isempty(i0)
                    x0 = x(i0);
                else
                    x0 = 0;
                end

                x(end+1) = x0 + optimValues.iteration;
                y(end+1) = optimValues.fval;

                set(obj.Line_c,'XData',x,'YData',y);
            %otherwise
            %    error('Unexpected method: "%s"',method);
        end
        drawnow

    case 'init'
        % Initialize the figure?
        drawnow
    otherwise
        fprintf('%s, %s\n',state,method);
end

end

% Error using matlab.graphics.axis.Axes/set
% The UIContextMenu's parent figure is not the same as the parent figure of this object.
% 
% Error in callAllOptimPlotFcns (line 128)
%         set(handle,'UIContextMenu', cmenu);
% 
% Error in barrier
% 
% Error in barrier
% 
% Error in barrier
% 
% Error in fmincon (line 862)
%     [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = barrier(funfcn,X,A,B,Aeq,Beq,l,u,confcn,options.HessFcn, ...
% 
% Error in approximateFixedKinematics (line 170)
%     [q,fval,exitflag,output] = fmincon(fcost_q,q,[],[],[],[],q_lb,q_ub,[],fOptions_q);
% 
% Error in SCRIPT_approximateFixedKinematics (line 29)
% [c,q,err] = approximateFixedKinematics(fcnH_e2o,H_e2o_ALL,options);
