function [p_w,err,p_w_i] = multiviewPointsToWorld(A_c2m,H_c2w,p_m)
% MULTIVIEWPOINTSTOWORLD recovers world-fixed points using point
% correspondeces between multiple cameras with known intrinsics and
% extrinsics.
%   p_w = multiviewPointsToWorld(p_m,A_c2m,H_c2w)
%
%   [p_w, err, p_w_i] = multiviewPointsToWorld(___)
%
%   Input(s)
%       A_c2m - N-element cell array containing intrinsics associated with
%               each camera
%           A_c2m{1} - (K+1)x(K+1) array defining intrinsics of Camera 1 
%           A_c2m{2} - (K+1)x(K+1) array defining intrinsics of Camera 2
%           ...
%           A_c2m{N} - (K+1)x(K+1) array defining intrinsics of Camera N
%       H_c2w - N-emelemnt cell array containing extrinsics referencing
%               each camera to a common world frame 
%           H_c2w{1} - (K+2)x(K+2) array defining extrinsics relating 
%                      Camera 1 to the world frame
%           H_c2w{2} - (K+2)x(K+2) array defining extrinsics relating 
%                      Camera 2 to the world frame
%           ...
%           H_c2w{N} - (K+2)x(K+2) array defining extrinsics relating 
%                      Camera N to the world frame
%       p_m   - N-element cell array containing point correspondences
%           p_m{1} - KxM array containing Camera 1 pixel coordinates for 
%                    M points 
%           p_m{2} - KxM array containing Camera 2 pixel coordinates for 
%                    M points
%           ...
%           p_m{N} - KxM array containing Camera N pixel coordinates for 
%                    M points 
%
%   Output(s)
%       p_w   - (K+1)xM array containing points relative to world frame
%       err   - NxM array containing euclidean distance between the set of
%               points p_w, and the points recovered by each individual 
%               camera
%           err(i,j) - distance between ith camera's recovery of point j,
%                      and p_w(:,j).
%       p_w_i - N-element cell array containing world-referenced point
%               coordinates recovered by each camera
%           p_w_i{1} - (K+1)xM array containing points from Camera 1 
%                      relative to world frame
%           p_w_i{2} - (K+1)xM array containing points from Camera 2 
%                      relative to world frame
%           ...
%           p_w_i{N} - (K+1)xM array containing points from Camera N 
%                      relative to world frame
%
%   Implementation Notes:
%       (1) This code is generically implemented to support ND cameras
%           where N >= 2 (a typical camera is 3D)
%       (2) This code does not fully check for valid intrinsics and 
%           extrinsics
%
%   Intrinsic, Extrinsic, and Point Example (3D)
%
%       Intrincis
%                      [ sx , sxy , u0 ]
%           A_c2m{i} = [  0 ,  sy , v0 ]
%                      [  0 ,   0 ,  1 ]
%
%       Extrinsics
%                      [ r11 , r12 , r13 , d1 ]
%           H_c2w{i} = [ r21 , r22 , r23 , d2 ]
%                      [ r31 , r32 , r33 , d3 ]
%                      [   0 ,   0 ,   0 ,  1 ]
%
%       Point
%           p_m{i}(:,j) = [x; y]
%
%   M. Kutzer, 27Sep2022, USNA

% TODO - allow cameraParameters to be used instead of intrinsics

% Updates:

%% Check input(s)
narginchk(3,3);

% Confirm all inputs are cell arrays
if ~iscell(A_c2m)
    error('Intrinsics must be specified as a cell array.');
end
if ~iscell(H_c2w)
    error('Extrinsics must be specified as a cell array.');
end
if ~iscell(p_m)
    error('Point correspondences must be specified as a cell array.');
end

% Confirm that all inputs are the same length
N = numel(A_c2m);
if numel(H_c2w) ~= N || numel(p_m) ~= N
    msg = sprintf([...
        'The number of intrinsic matrices, extrinsic matrices, and ',...
        'sets of points must match\n',...
        '\tnumel(A_c2m) = %d\n',...
        '\tnumel(H_c2m) = %d\n',...
        '\tnumel(p_m)   = %d\n',...
        'These values must all be the same.'],...
        numel(A_c2m),numel(H_c2w),numel(p_m) );
    error(msg);
end

% Check points
[K,M,tmp] = size(p_m{1});
if tmp ~= 1
    error('p_m{i} must be a KxM (2D array only).');
end

if any( cellfun(@(x)size(x,1) ~= K || size(x,2) ~= M || size(x,3) ~= 1,p_m) )
    error('Each element of p_m must be a KxM array.');
end

% Check dimensions of intrinsics
if any( cellfun(@(x)size(x,1) ~= (K+1) || size(x,2) ~= (K+1) || size(x,3) ~= 1,A_c2m) )
    error('Each element of the provided intrinsics must be a %dx%d array for the points provided.',(K+1),(K+1));
end

% Check dimensions of extrinsics
if any( cellfun(@(x)size(x,1) ~= (K+2) || size(x,2) ~= (K+2) || size(x,3) ~= 1,H_c2w) )
    error('Each element of the provided extrinsics must be a %dx%d array for the points provided.',(K+2),(K+2));
end

% Check last row of intrinsics
lastRow = zeros(1,K+1);
lastRow(end) = 1;
if any( cellfun(@(x)~isequal(x(end,:),lastRow),A_c2m) )
    msg = sprintf([...
        'Each element of the provided intrinsics must have a last row ',...
        'equal to [']);
    numMsg = sprintf(' %d ',lastRow);
    msg = sprintf('%s%s].',msg,numMsg);
    error(msg);
end

% Check last row of extrinsics
lastRow = zeros(1,K+2);
lastRow(end) = 1;
if any( cellfun(@(x)~isequal(x(end,:),lastRow),H_c2w) )
    msg = sprintf([...
        'Each element of the provided extrinsics must have a last row ',...
        'equal to [']);
    numMsg = sprintf(' %d ',lastRow);
    msg = sprintf('%s%s].',msg,numMsg);
    error(msg);
end

%% Calculate scaled coordinates relative to the camera frame
tilde_p_c = repmat({nan(K+1,M)},1,N);
for i = 1:N
    % TODO - allow cameraParameters to be used instead of intrinsics

    % Calculate inverse of intrinsics
    A_m2c = A_c2m{i}^(-1);

    % Check for finite inverse
    if ~all( isfinite(A_m2c) )
        error('A_c2m{%d} is not invertible',i);
    end
    
    % Make pixel coordinates homogeneous
    p_m{i}(end+1,:) = 1;

    % Calculate unscaled points relative to camera frame
    tilde_p_c{i} = A_m2c * p_m{i};
end

%% Calculate matrices to solve multiview problem
% Z = pinv(B)*D

% Calculate H_ci2cj & D
%   [i,j] = [2 ,  1 ]
%   [i,j] = [3 ,  2 ]
%   [i,j] = [4 ,  3 ]
%   ...
%   [i,j] = [N , N-1]
%   [i,j] = [1 ,  N ]
D = nan((K+1)*N,1);
H_ci2cj = repmat({nan(K+2,K+2)},1,N);
for i = 2:N
    % Define (i-1)th relative transform
    H_ci2cj{i-1} = invSE(H_c2w{i-1},true)*H_c2w{i};

    % Populate (i-1)th block of D
    idx = (i-2)*(K+1)+(1:(K+1));
    D(idx,1) = H_ci2cj{i-1}( 1:(end-1),end );
end
H_ci2cj{N} = invSE(H_c2w{N},true)*H_c2w{1};
% Populate Nth block of D
idx = (N-1)*(K+1)+(1:(K+1));
D(idx,1) = H_ci2cj{N}( 1:(end-1),end );

% Calculate B
B = repmat({zeros((K+1)*N,N)},1,M);
for k = 1:M
    for i = 1:N
        idx = (i-1)*(K+1)+(1:(K+1));
        for j = 1:N

            if i == j
                B{k}(idx,j) = tilde_p_c{j}(:,k);
                continue
            end

            if i+1 == j
                B{k}(idx,j) = -H_ci2cj{i}(1:(end-1),1:(end-1))*tilde_p_c{j}(:,k);
                continue
            end

            if i == N && j == 1
                B{k}(idx,j) = -H_ci2cj{i}(1:(end-1),1:(end-1))*tilde_p_c{j}(:,k);
                continue
            end

        end
    end
end

% Calculate depth
Z = nan(N,M);
for k = 1:M
    Z(:,k) = pinv(B{k})*D;
end

%% Calculate scaled points & points relative to world frame
p_c = nan(K+2,M);
p_w_i = repmat({nan(K+2,M)},1,N);
for i = 1:N
    p_c(1:(K+1),:) = Z(i,:).*tilde_p_c{i}; % Points relative to camera i
    p_c(K+2,:) = 1;                        % Make points homogeneous
    p_w_i{i} = H_c2w{i}*p_c;               % Reference points to world frame
end

%% Combine points using mean
p_w = mean(reshape(cell2mat(p_w_i),M*(K+2),N), 2);
p_w = reshape(p_w,K+2,M);
p_w(end,:) = [];    % Remove row of ones

%% Calculate errors
if nargout > 1
    err = nan(N,M);
    for i = 1:N
        p_w_i{i}(end,:) = []; % Remove row of ones

        % Calculate difference
        dp_w = p_w_i{i}(1:(end-1),:) - p_w;
        % Calculate Euclidean distance
        err(i,:) = sqrt( sum(dp_w.^2,1) );
    end
end