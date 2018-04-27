function funcJ = calculateBodyFixeddJacobian(q,H,varargin)
%calculateBodyFixeddJacobian calculates the manipulator Jacobian relative 
%to body-fixed end-effector frame
%   calculateBodyFixeddJacobian(q,H) calculates the manipulator Jacobian  
%   associated with forward kinematics defined by transformation "H", and  
%   joint variables "q". Note, H can be an element of SE(2) or SE(3), and H 
%   and q must be symbolic.
%   This function returns an anonymous function, "funcJ" with input 
%   vector q.
%
%   calculateBodyFixeddJacobian(q,H,'file','filename') calculates the  
%   manipulator Jacobian associated with forward kinematics defined by  
%   transformation "H", and joint variables "q". Note, H can be an element 
%   of SE(2) or SE(3), and H and q must be symbolic.
%   This function saves a function as "filename" in the current directory,
%   with input vector q. This function returns the associated function
%   handle.
%
%   M. Kutzer 24Dec2014, USNA

%% Check inputs 
if nargin < 2
    error('Both "H" and "q" must be specified.')
end
if ~strcmpi( class(H), 'sym')  || ~strcmpi( class(q), 'sym')
    error('"H" and "q" must be symbolic variables.');
end
%TODO - check for properties of SE(2) and SE(3)
if size(H,1) == 4 && size(H,2) == 4 && ismatrix(H)
    dim = 3;
end
if size(H,1) == 3 && size(H,2) == 3 && ismatrix(H)
    dim = 2;
end

%% Check for custom functions
%TODO - check for vee.m

%% Calculate translation Jacobian
X = H(1:dim,dim+1); % translation associated with forward kinematics

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
R = H(1:dim,1:dim);     % rotation associated with forward kinematics
JT = transpose(R)*JT;   % transform velocity vector to body fixed frame

delete(h);
drawnow;

%% Calculate rotation Jacobian
%R = H(1:dim,1:dim); % rotation associated with forward kinematics

h = waitbar(0,'Calculating rotation portion of Jacobian...');
fprintf('Calculating rotation portion of Jacobian...\n');
%n = numel(q);
iter = 0;
for j = 1:n
    fprintf('Differentiating with respect to Joint %d...',j);
    dR = diff(R,q(j));
    fprintf('DONE\n');
    fprintf('Calculating se(n)...');
    %M = transpose(R)*dR;   % Jacobian relative to base
    M = transpose(dR)*R;    % Jacobian relative to end-effector? %CHECK ME
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