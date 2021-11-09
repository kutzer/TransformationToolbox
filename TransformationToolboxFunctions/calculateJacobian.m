function funcJ = calculateJacobian(q,H_e2o,varargin)
%calculateJacobian calculates the manipulator Jacobian relative to the base
%frame of the manipulator. 
%   funcJ = calculateJacobian(q,H_e2o) calculates the manipulator Jacobian  
%   associated with forward kinematics defined by transformation "H_e2o",   
%   and joint variables "q". Note, H_e2o can be an element of SE(2) or  
%   SE(3), and H_e2o and q must be symbolic. This function returns an 
%   anonymous function, "funcJ" with input vector q.
%
%   funcJ = calculateJacobian(___,'file','filename') this function saves a 
%   function as "filename" in the current directory, with input vector q. 
%   This function returns the associated function handle.
%
%   Input(s)
%       q     - Nx1 array containing N symbolic joint variables
%       H_e2o - 3x3 element of SE(2) or 4x4 element of SE(3) defining the
%               forward kinematics of the manipulator using symbolic terms
%               for joint variables (all must be contained in q)
%
%   Output(s)
%       funcJ - function handle for the Jacobian such that J = funcJ(q).
%               The returned J matrix will be a 6xN array.
%
%   Example(s)
%       % Given H_e2o_now, q_now, and H_e2o_des; find delta_q
%       
%       %   This assumes H_e2o_now and H_e2o_des are *close together*
%       delta_H_e2o = invSE(H_e2o_now) * H_e2o_des;
%       % Calculate delta_X
%       delta_X(1:3,1) = delta_H_e2o(1:3,4);
%       delta_X(4:6,1) = vee( logm(delta_H_e2o(1:3,1:3)),'fast');
%       % Calculate delta_q
%       delta_q = pinv( funcJ(q_now) )*delta_X;
%       % Calculate updated q
%       q = q_now + delta_q;
%
%   M. Kutzer 10Oct2014, USNA

% Updates:
%   09Nov2021 - Updated documentation
%% Check inputs 
if nargin < 2
    error('Both "H" and "q" must be specified.')
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

%% Check for custom functions
%TODO - check for vee.m

%% Calculate translation Jacobian
X = H_e2o(1:dim,dim+1); % translation associated with forward kinematics

h = waitbar(0,'Calculating translation portion of Jacobian...');
fprintf('Calculating translation portion of Jacobian...\n');
m = numel(X);
n = numel(q);
iter = 0;
for i = 1:m
    fprintf('Calculating for dimension %d...\n',i);
    for j = 1:n
        fprintf('Differentiating with respect to Joint %d...',j);
        JT(i,j) = simplify( diff(X(i),q(j)) );
        iter = iter+1;
        waitbar(iter/(m*n),h)
        fprintf('DONE\n');
    end
end
delete(h);
drawnow;

%% Calculate rotation Jacobian
R = H_e2o(1:dim,1:dim); % rotation associated with forward kinematics

h = waitbar(0,'Calculating rotation portion of Jacobian...');
fprintf('Calculating rotation portion of Jacobian...\n');
%n = numel(q);
iter = 0;
for j = 1:n
    fprintf('Differentiating with respect to Joint %d...',j);
    dR = diff(R,q(j));
    fprintf('DONE\n');
    fprintf('Calculating se(n)...');
    M = transpose(R)*dR;
    fprintf('DONE\n');
    fprintf('Vectorizing se(n)...');
    v = vee(M,'fast');
    fprintf('DONE\n');
    fprintf('Combining result...');
    JR(:,j) = v;
    iter = iter+1;
    waitbar(iter/n,h)
    fprintf('DONE\n');
end
delete(h);
drawnow;

%% Deal with constant syms
%TODO compensate for constant syms (e.g. l1, l2, l3)
%% Combine to create full Jacobian using an anonymous function
if isempty(varargin)
    fprintf('Creating function handle...');
    funcJ = matlabFunction([JT; JR],'vars',{q});
    fprintf('DONE\n');
    return
end

%% Try to create a file
if ~isempty(varargin)
    fprintf('Saving function "%s.m" to Current Directory...',varargin{2});
    funcJ = matlabFunction([JT; JR],varargin{1},varargin{2},'vars',{q});
    fprintf('SAVED\n');
    return
end