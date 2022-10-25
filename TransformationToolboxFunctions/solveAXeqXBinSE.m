function [X,A,B] = solveAXeqXBinSE(A,B,varargin)
% SOLVEAXEQXBINSE solves the AX = XB equation assuming A, B, and X \in SE(N) 
% given two or more pairs of A and B following [1]
%   X = SOLVEAXEQXBINSE(A,B) provides the AX = XB least squares solution 
%   given two or more pairs of A and B. 
%   
%   X = SOLVEAXEQXBINSE(___,ZERO)
%
%   X = SOLVEAXEQXBINSE(___,fast)
%
%   [X,A_used,B_used] = solveAXeqXBinSE(___)
%
%   Input(s)
%       A - k-element cell array containing values of SE(N)
%       B - k-element cell array containing values of SE(N)
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not   
%              specified, a default value is used.
%       fast - [OPTIONAL] true/false logical value indicating whether to
%              skip checking SE(N).
%                fast = true    - Skip checking if H \in SE(N)
%                fast = [false] - Check if H \in SE(N)
%
%   Output(s)
%       X      - (N+1)x(N+1) matrix containing the least squares solution  
%                to AX = XB
%       A_used - m-element cell array containing "corrected" values of A
%               (see validCorrespondenceSE.m)
%       B_used - m-element cell array containing "corrected" values of B
%               (see validCorrespondenceSE.m)
%
%   References:
%       [1] Park, Frank C., and Bryan J. Martin. "Robot sensor calibration:
%           solving AX= XB on the Euclidean group." IEEE Transactions on 
%           Robotics and Automation 10.5 (1994): 717-721.
%
%   M. Kutzer, 09Apr2021, USNA

% Update(s)
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast
%   25Oct2022 - Preallocate C/d and show waitbar for large data sets
%% Default options
ZERO = [];
fast = false;

%% Check input(s)
narginchk(2,4);
if ~iscell(A) || ~iscell(B)
    error('A and B must be specified as k-element cell arrays.');
end

if numel(A) ~= numel(B)
    error('A and B must be specified in pairs (i.e. numel(A) == numel(B)).');
end

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

%% Correct A/B pairs
str = 'AB';
if ~fast
    [A,B,info] = validCorrespondenceSE(A,B,ZERO);
    
    % Display removed pair information and altered transforms
    for j = 1:size(info.RemoveBin,2)
        % Display removed pair information
        if info.RemoveBin(1,j)
            fprintf('REMOVED PAIR: A{%d} - %s, B{%d} - %s\n',...
                j,info.RemoveMsg{1,j},j,info.RemoveMsg{2,j});
            continue
        end
        
        % Display altered transform information
        for i = 1:size(info.RemoveBin,1)
            if ~isempty(info.RemoveMsg{i,j})
                fprintf('ALTERED TRANSFORM: %s{%d} - %s\n',...
                    str(i),j,info.RemoveMsg{i,j});
            end
        end
    end

end

%% Calculate M to solve for R_X

% Define dimensions
n = size(A{1},1);

% Initialize M
M = zeros(n-1,n-1);

fastLocal = true;
k = numel(A);

% Provide status for large sets of data
if k >= 10000
    tfStatus = true;
    nStatus = 50;
    iStatus = unique( round(linspace(1,k,nStatus)) );
    msg = 'Solving AX=XB rotation...';
    wb = waitbar(0,msg,'Name','solveAXeqXB.m');
    dt_Status = now;
    t0_Status = tic;
else
    tfStatus = false;
end

% Solve rotation
for i = 1:k
    % Calculate vectorized forms of so
    R_Ai = A{i}(1:(n-1),1:(n-1));    % ith rotation for A
    R_Bi = B{i}(1:(n-1),1:(n-1));    % ith rotation for B

    % TODO - check Tr(R_A) ~= -1 & Tr(R_B) ~= -1
    a = vee( logSO(R_Ai,fastLocal), fastLocal );  % vectorized form of so
    b = vee( logSO(R_Bi,fastLocal), fastLocal );  % vectorized form of so
    
    % Update M (page 720)
    M = M + b*a.';

    % Update status
    if tfStatus
        if any(i == iStatus)
            prc = i/k;
            ti_Status = toc(t0_Status);
            tf_Status = (1-prc)*ti_Status/prc;
            msg_i = sprintf('~ %d seconds remaining',round(tf_Status));
            waitbar(prc,wb,{msg,msg_i});
            drawnow;
        end
    end
end

% Solve for rotation (page 720)
R_X = ( ( (M.')*M )^(-1/2) )*(M.');

% Delete status waitbar
if tfStatus
    delete(wb);
    drawnow
end

%% Calculate C and d to solve for translation b_X

if tfStatus
    msg = 'Solving AX=XB translation...';
    wb = waitbar(0,msg,'Name','solveAXeqXB.m');
    dt_Status = now;
    t0_Status = tic;
end

%C = [];
%d = [];
C = zeros(k*(n-1),(n-1));
d = zeros(k*(n-1),1);
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
    %C = [C; Ci];
    %d = [d; di];
    C( ( (i-1)*(n-1)+(1:(n-1)) ),1:(n-1) ) = Ci;
    d( ( (i-1)*(n-1)+(1:(n-1)) ),1       ) = di;

    % Update status
    if tfStatus
        if any(i == iStatus)
            prc = i/k;
            ti_Status = toc(t0_Status);
            tf_Status = (1-prc)*ti_Status/prc;
            msg_i = sprintf('~ %d seconds remaining',round(tf_Status));
            waitbar(prc,wb,{msg,msg_i});
            drawnow;
        end
    end
end
% Solve for b_X (page 720)
b_X = ( ((C.')*C)^(-1) )*(C.')*d;

% Delete status waitbar
if tfStatus
    delete(wb);
    drawnow;
end
%% Compile output
X = eye(n);
X(1:(n-1),1:(n-1)) = R_X;
X(1:(n-1),n)       = b_X;