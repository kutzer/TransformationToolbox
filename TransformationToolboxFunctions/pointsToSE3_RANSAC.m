function [H_q2p,inlierIdx,errInlier,errAll] = ...
    pointsToSE3_RANSAC(X_q,X_p,sampleSize,maxDistance,varargin)
% POINTSTOSE3_RANSAC finds the best fit rigid body between two sets of
% point correspondences using random sample consensus (RANSAC).
%   [H_q2p,inlierIdx,errInlier,errAll] = ...
%       pointsToSE3_RANSAC(X_q,X_p,sampleSize,maxDistance)
%
%   ___ = pointsToSE3_RANSAC(___, Name, Value) specifies additional 
%   name-value pair arguments described from ransac.m
%
%   Input(s)
%           X_q - 3xN array of points referenced to frame q
%           X_p - 3xN array of points referenced to frame p
%    sampleSize - scalar defining minimum sample size
%   maxDistance - scalar defining maximum distance of inliers
%
%   Output(s)
%       H_q2p - rigid body transformation defining frame q relative to
%               frame p
%   inlierIdx - 
%
%   See also pointsToSE3 ransac
%
%   M. Kutzer, 05Dec2023, USNA

%% Check input(s)
narginchk(4,inf);

N = size(X_p,2);
if size(X_q,2) ~= N
    error('p and q must be the same dimension.');
end

if size(X_q,1) ~= 3
    error('The first input must be a 3xN numeric array.');
end

if size(X_p,1) ~= 3
    error('The second input must be a 3xN numeric array.');
end

% TODO - check input(s)

%% Solve RANSAC
data = packData(X_q,X_p);

fitFcn = @(data) fitTformFcn(data);
distFcn = @(dataFit,data) evalTformFcn(dataFit,data);

[H_q2p,inlierIdx] = ransac(data,fitFcn,distFcn,sampleSize,maxDistance,...
    varargin{:});

%% Package output(s)
idxInlier = inlierIdx.';
errInlier = distFcn(H_q2p,data(inlierIdx,:));
errAll = distFcn(H_q2p,data);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Pack data
function data = packData(X_q,X_p)

data = [X_q; X_p].';

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Unpack data
function [X_q,X_p] = unpackData(data)

X_q = data(:,1:3).';
X_p = data(:,4:6).';

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Fit transformation to data
function H_q2p = fitTformFcn(data)

[X_q,X_p] = unpackData(data);
H_q2p = pointsToSE3(X_q,X_p);

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Evaluate tform using data
function d = evalTformFcn(H_q2p,data)

[X_q,X_p] = unpackData(data);
sX_p = H_q2p*[X_q; ones(1,size(X_q,2))];

dX_p = sX_p(1:3,:) - X_p(1:3,:);

d = sum( dX_p.^2 , 1 ).';

end