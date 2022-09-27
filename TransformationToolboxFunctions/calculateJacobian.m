function funcJ = calculateJacobian(q,H_e2o,varargin)
%CALCULATEJACOBIAN calculates the world or body-fixed manipulator Jacobian
%   funcJ = calculateJacobian(q,H_e2o)
%
%   funcJ = calculateJacobian(q,H_e2o,showStatus)
%
%   funcJ = calculateJacobian(___,'Name','Value')
%
%   Input(s)
%       q     - Nx1 array containing N symbolic joint variables
%       H_e2o - 3x3 element of SE(2) or 4x4 element of SE(3) defining the
%               forward kinematics of the manipulator using symbolic terms
%               for joint variables (all must be contained in q)
%       showStatus - [Default, true] shows the status in the command prompt
%                    if true, hides status otherwise. Note that saved
%                    function status is shown regardless.
%
%       Name-Value Arguments
%           Specify optional comma-separated pairs of Name,Value arguments.
%           Name is the argument name and Value is the corresponding value.
%
%              Name | Value [default]
%       'Reference' | { ['Base'], 'End-effector' }
%           'Order' | { ['TranslationRotation'], 'RotationTranslation' }
%       'Constants' | Mx1 symbolic array containing kinematic constants
%        'Filename' | character array specifying Jacobian function name
%
%       NOTE: 'world' and 'body' values for 'reference' are also acceptable
%
%   Output(s)
%       funcJ - function handle for the Jacobian such that J = funcJ(q)
%               returns a 6xN array defined:
%
%       Order:
%           'TranslationRotation'
%               funcJ(q) = [ Jtran(q) ]
%                          [ Jrot(q)  ]
%           'RotationTranslation'
%               funcJ(q) = [ Jrot(q)  ]
%                          [ Jtran(q) ]
%
%   M. Kutzer 10Oct2014, USNA

% TODO - Add coupled/decoupled Jacobian option! 

% Updates:
%   09Nov2021 - Updated documentation
%   23Feb2022 - Revised to correctly differentiate between body/world
%               referenced Jacobians
%   24Feb2022 - Revised to allow user specification of translation/rotation
%               vs rotation/translation stacking of Jacobian
%   28Feb2022 - Documentation update and added comments for saved function
%   03Mar2022 - Added showStatus input option
%   12Apr2022 - Added support of symbolic constant kinematic terms
%   27Apr2022 - Added symbolic constant kinematic terms to saved function
%               option
%% Check inputs
if nargin < 2
    error('Both "H_e2o" and "q" must be specified.')
end
if ~strcmpi( class(H_e2o), 'sym')  || ~strcmpi( class(q), 'sym')
    error('"H" and "q" must be symbolic variables.');
end
%TODO - check for properties of SE(2) and SE(3)
if size(H_e2o,1) == 4 && size(H_e2o,2) == 4 && ismatrix(H_e2o)
    dim = 3;
end
if size(H_e2o,1) == 3 && size(H_e2o,2) == 3 && ismatrix(H_e2o)
    dim = 2;
end

%% Parse name/value pairs
% Define default values
c = [];
filename = [];
reference = 'base';
order = 'translationrotation';
if nargin > 2
    nV = numel(varargin);
    if ~ischar(varargin{1})
        showStatus = varargin{1};
        varargin(1) = [];
        nV = numel(varargin);
    else
        showStatus = false;
    end
    
    if nV/2 ~= round( nV/2 )
        error('Name-Value Arguments must be specified in pairs.')
    end
    
    for i = 1:2:nV
        Name = varargin{i};
        Value = varargin{i+1};
        switch lower(Name)
            case 'reference'
                reference = lower( Value );
            case 'order'
                order = lower( Value );
            case 'constants'
                c = Value;
                % TODO - check if values are constant!
            case 'filename'
                filename = Value;
            otherwise
                error('"%s" is not a valid property name.')
        end
    end
end

% Check reference
switch reference
    case 'base'
        % Jacobian is world-referenced
        isWorld = true;
    case 'world'
        % Jacobian is world-referenced
        isWorld = true;
    case 'end-effector'
        % Jacobian is body-fixed
        isWorld = false;
    case 'body'
        % Jacobian is body-fixed
        isWorld = false;
    otherwise
        error('Value "%s" is not defined for the property name "Reference".',reference);
end

% Check order
switch order
    case 'translationrotation'
        isTransRot = true;
    case 'rotationtranslation'
        isTransRot = false;
    case 'tr'
        isTransRot = true;
    case 'rt'
        isTransRot = false;
    otherwise
        error('Value "%s" is not defined for the property name "Order".',order);
end

% TODO - check filename

%% Check for custom functions
%TODO - check for vee.m

%% Parse rotation & translation from rigid body transform
R_e2o = H_e2o(1:dim,1:dim); % rotation associated with forward kinematics
d_e2o = H_e2o(1:dim,dim+1); % translation associated with forward kinematics

%% Calculate translation Jacobian
if showStatus
    h = waitbar(0,'Calculating translation portion of Jacobian...');
    fprintf('Calculating translation portion of Jacobian...\n');
end
%m = numel(d_e2o);
n = numel(q);
iter = 0;
% for i = 1:m
%     fprintf('\tCalculating for dimension %d...\n',i);
%     for j = 1:n
%         fprintf('\t\tDifferentiating with respect to Joint %d...',j);
%         Jtran(i,j) = simplify( diff(d_e2o(i),q(j)) );
%         iter = iter+1;
%         waitbar(iter/(m*n),h)
%         fprintf('DONE\n');
%     end
% end
for j = 1:n
    if showStatus
        fprintf('\tDifferentiating with respect to Joint %d...',j);
    end
    Jtran(:,j) = simplify( diff(d_e2o,q(j)) );
    if showStatus
        iter = iter+1;
        waitbar(iter/n,h)
        fprintf('DONE\n');
    end
end
if showStatus
    delete(h);
    drawnow;
end

% Update for world/body referenced
if ~isWorld
    Jtran = transpose(R_e2o)*Jtran;
end

%% Calculate rotation Jacobian
if showStatus
    h = waitbar(0,'Calculating rotation portion of Jacobian...');
    fprintf('Calculating rotation portion of Jacobian...\n');
end
%n = numel(q);
iter = 0;
for j = 1:n
    if showStatus
        fprintf('\tDifferentiating with respect to Joint %d...',j);
    end
    dR = diff(R_e2o,q(j));
    if showStatus
        fprintf('DONE\n');
        fprintf('\t\tCalculating so(n)...');
    end
    if isWorld
        K = dR*transpose(R_e2o);
    else
        K = transpose(R_e2o)*dR;
    end
    if showStatus
        fprintf('DONE\n');
        fprintf('\t\tVectorizing so(n)...');
    end
    k = vee(K,'fast');
    if showStatus
        fprintf('DONE\n');
        fprintf('\t\tCombining result...');
    end
    Jrot(:,j) = k;
    if showStatus
        iter = iter+1;
        waitbar(iter/n,h)
        fprintf('DONE\n');
    end
end
if showStatus
    delete(h);
    drawnow;
end

%% Deal with constant syms
%TODO compensate for constant syms (e.g. l1, l2, l3)

%% Combine Jacobian
if isTransRot
    J = [Jtran; Jrot];
else
    J = [Jrot; Jtran];
end

%% Combine to create full Jacobian using an anonymous function
if showStatus
    fprintf('Creating function handle...');
end
if ~isempty(c)
    funcJ = matlabFunction(J,'vars',{q,c});
else
    funcJ = matlabFunction(J,'vars',{q});
end
if showStatus
    fprintf('DONE\n');
end

%% Create a saved MATLAB function if a filename is supplied
if ~isempty(filename)
    % Define function comments
    if ~isempty(c)
        fcnIn = 'q, c';
    else
        fcnIn = 'q';
    end

    if isWorld
        ref = 'world';
    else
        ref = 'body';
    end
    
    if isTransRot
        jStack = sprintf('J(%s) = [ Jtran(%s); Jrot(%s) ]',fcnIn,fcnIn,fcnIn);
    else
        jStack = sprintf('J(%s) = [ Jrot(%s); Jtran(%s) ]',fcnIn,fcnIn,fcnIn);
    end
    
    str = sprintf([...
        '%s calculates a %s-referenced manipulator Jacobian\n'],...
        upper(filename),ref);
    str = sprintf(['%s',...
        '\tJ = %s(%s)\n\n'],...
        str,filename,fcnIn);
    str = sprintf(['%s',...
        '\tInput(s)\n',...
        '\t\tq - %dx1 array defining joint configuration\n'],...
        str,numel(q));
    if ~isempty(c)
        str = sprintf(['%s',...
        '\t\tc - %dx1 array defining fixed parameters\n'],...
        str,numel(c));
    end
    str = sprintf('%s\n',str);
    str = sprintf(['%s',...
        '\tOutput(s)\n',...
        '\t\tJ - %dx%d array defining %s-referenced Jacobian for \n',...
        '\t\t    joint configuration q\n\n',...
        '\t\t    %s\n\n'],...
        str,size(J,1),size(J,2),ref,jStack);
    str = sprintf(['%s',...
        '\tFunction created using "caculateJacobian.m"\n',...
        '\t%s\n'],...
        str,datestr(now));
    
    fprintf('\n-------- DETAILED FUNCTION DOCUMENTATION (copy/paste) --------\n')
    fprintf('%s',str);
    fprintf('--------------------------------------------------------------\n\n')
    
    % Save function
    fprintf('Saving function "%s.m" to Current Directory...',filename);
    if ~isempty(c)
        funcJ = matlabFunction(J,'file',filename,'vars',{q,c},...
            'Comments',sprintf('calculateJacobian: %s-referenced, %s',ref,jStack));
    else
        funcJ = matlabFunction(J,'file',filename,'vars',{q},...
            'Comments',sprintf('calculateJacobian: %s-referenced, %s',ref,jStack));
    end
    % TODO - updated to include detailed comments if/when matlabFunction.m
    %        improves support
    %funcJ = matlabFunction(J,'file',filename,'vars',{q},...
    %    'Comments',str);
    fprintf('SAVED\n');
    return
end