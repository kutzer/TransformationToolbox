function X = solveAXeqXBinSE(A,B)
% SOLVEAXEQXBINSE solves the AX = XB equation assuming A, B, and X \in SE(N) 
% given two or more pairs of A and B following [1]
%   X = SOLVEAXEQXBINSE(A,B) provides the AX = XB least squares solution 
%   given two or more pairs of A and B. 
%
%   Input(s)
%       A - k-element cell array containing values of SE(N)
%       B - k-element cell array containing values of SE(N)
%
%   Output(s)
%       X - (N+1)x(N+1) matrix containing the least squares solution to 
%           AX = XB
%
%   References:
%       [1] Park, Frank C., and Bryan J. Martin. "Robot sensor calibration:
%           solving AX= XB on the Euclidean group." IEEE Transactions on 
%           Robotics and Automation 10.5 (1994): 717-721.
%
%   M. Kutzer, 09Apr2021, USNA

%% Check input(s)
narginchk(2,2);
if ~iscell(A) || ~iscell(B)
    error('A and B must be specified as k-element cell arrays.');
end

if numel(A) ~= numel(B)
    error('A and B must be specified in pairs (i.e. numel(A) == numel(B)).');
end

%% Calculate M to solve for R_X
[n,m,o] = size(A{1});
if o ~= 1
    error('Cell elements must be NxNx1');
end

% Initialize M
M = zeros(n-1,n-1);

k = numel(A);
for i = 1:k
    % Check elements
    % -> Check matrix dimensions
    [niA,miA] = size(A{i});
    [niB,miB] = size(B{i});
    if niA ~= miA
        error('A(%d) is not square.',i);
    end
    if niB ~= miB
        error('B(%d) is not square.',i);
    end
    if (niA ~= n) || (miA ~= m) || (niB ~= n) || (miB ~= m)
        error('All matrices must be the same dimension');
    end
    % -> Valid elements of SE
    [tfA,msgA] = isSE(A{i});
    [tfB,msgB] = isSE(B{i});
    if ~tfA
        fprintf('A(%d) does not appear to be an element of SE(%d):\n',...
            i,n-1);
        fprintf('\t%s\n',msgA);
    end
    if ~tfB
        fprintf('B(%d) does not appear to be an element of SE(%d):\n',...
            i,n-1);
        fprintf('\t%s\n',msgB);
    end
    
    % Calculate vectorized forms of so
    R_Ai = A{i}(1:(n-1),1:(n-1));    % ith rotation for A
    R_Bi = B{i}(1:(n-1),1:(n-1));    % ith rotation for B
    % TODO - check Tr(R_A) ~= -1 & Tr(R_B) ~= -1
    a = vee(logm(R_Ai),'fast');  % vectorized form of so
    b = vee(logm(R_Bi),'fast');  % vectorized form of so
    
    % Update M (page 720)
    M = M + b*a.';
end
% Solve for rotation (page 720)
R_X = ( ( (M.')*M )^(-1/2) )*(M.');

%% Calculate C and d to solve for translation b_X
C = [];
d = [];
for i = 1:k
    % Isolate Rotation of A
    R_Ai = A{i}(1:(n-1),1:(n-1));
    % Update C (page 720)
    Ci = eye(n-1) - R_Ai;
    
    % Isolate translations
    b_Ai = A{i}(1:(n-1),n);
    b_Bi = B{i}(1:(n-1),n);
    % Update d (page 720)
    di = b_Ai - R_X*b_Bi;
    
    % Compile C and d
    C = [C; Ci];
    d = [d; di];
end
% Solve for b_X (page 720)
b_X = ( ((C.')*C)^(-1) )*(C.')*d;

%% Compile output
X = eye(n);
X(1:(n-1),1:(n-1)) = R_X;
X(1:(n-1),n) = b_X;