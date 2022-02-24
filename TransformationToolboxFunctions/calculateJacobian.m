function funcJ = calculateJacobian(q,H_e2o,varargin)
%CALCULATEJACOBIAN calculates the world or body-fixed manipulator Jacobian
%   funcJ = calculateJacobian(q,H_e2o)
%
%   funcJ = calculateJacobian(___,'Name','Value')
%
%   Input(s)
%       q     - Nx1 array containing N symbolic joint variables
%       H_e2o - 3x3 element of SE(2) or 4x4 element of SE(3) defining the
%               forward kinematics of the manipulator using symbolic terms
%               for joint variables (all must be contained in q)
%
%       Name-Value Arguments
%           Specify optional comma-separated pairs of Name,Value arguments. 
%           Name is the argument name and Value is the corresponding value.
%
%              Name | Value [default]
%       'Reference' | { ['Base'], 'End-effector' }
%           'Order' | { ['TranslationRotation'], 'RotationTranslation' }
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
%               funcJ(q) = [ Jrot(q)  ]
%                          [ Jtran(q) ]
%           'RotationTranslation'
%               funcJ(q) = [ Jrot(q)  ]
%                          [ Jtran(q) ]
%
%   M. Kutzer 10Oct2014, USNA

% Updates:
%   09Nov2021 - Updated documentation
%   23Feb2022 - Revised to correctly differentiate between body/world
%               referenced Jacobians
%   24Feb2022 - Revised to allow user specification of translation/rotation
%               vs rotation/translation stacking of Jacobian

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
filename = [];
reference = 'base';
order = 'translationrotation';
if nargin > 2
    nV = numel(varargin);
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
h = waitbar(0,'Calculating translation portion of Jacobian...');
fprintf('Calculating translation portion of Jacobian...\n');
m = numel(d_e2o);
n = numel(q);
iter = 0;
for i = 1:m
    fprintf('\tCalculating for dimension %d...\n',i);
    for j = 1:n
        fprintf('\t\tDifferentiating with respect to Joint %d...',j);
        Jtran(i,j) = simplify( diff(d_e2o(i),q(j)) );
        iter = iter+1;
        waitbar(iter/(m*n),h)
        fprintf('DONE\n');
    end
end
delete(h);
drawnow;

% Update for world/body referenced
if ~isWorld
    Jtran = transpose(R_e2o)*Jtran;
end

%% Calculate rotation Jacobian
h = waitbar(0,'Calculating rotation portion of Jacobian...');
fprintf('Calculating rotation portion of Jacobian...\n');
%n = numel(q);
iter = 0;
for j = 1:n
    fprintf('\tDifferentiating with respect to Joint %d...',j);
    dR = diff(R_e2o,q(j));
    fprintf('DONE\n');
    fprintf('\tCalculating so(n)...');
    if isWorld
        K = dR*transpose(R_e2o);
    else
        K = transpose(R_e2o)*dR;
    end
    fprintf('DONE\n');
    fprintf('\tVectorizing so(n)...');
    k = vee(K,'fast');
    fprintf('DONE\n');
    fprintf('\tCombining result...');
    Jrot(:,j) = k;
    iter = iter+1;
    waitbar(iter/n,h)
    fprintf('DONE\n');
end
delete(h);
drawnow;

%% Deal with constant syms
%TODO compensate for constant syms (e.g. l1, l2, l3)

%% Combine Jacobian
if isTransRot
    J = [Jtran; Jrot];
else
    J = [Jrot; Jtran];
end

%% Combine to create full Jacobian using an anonymous function
fprintf('Creating function handle...');
funcJ = matlabFunction(J,'vars',{q});
fprintf('DONE\n');

%% Try to create a file
if ~isempty(filename)
    fprintf('Saving function "%s.m" to Current Directory...',filename);
    funcJ = matlabFunction(J,'file',filename,'vars',{q});
    fprintf('SAVED\n');
    return
end