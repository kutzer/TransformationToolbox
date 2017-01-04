function [idx,err] = findPointCorrespondence(q,p)
% FINDPOINTCORRESPONDENCE find the indeces associated with the best fit
% correspondence between two sets of points.
%   idx = FINDPOINTCORRESPONDENCE(q,p) returns the linear indices
%   corresponding to the least squares best fit correspondence between q
%   and p assuming they are rigidly defined.
%
%   Given two sets of N corresponding points:
%       p = [p_1,p_2,p_3...p_N], p_i - 3x1 
%       q = [q_1,q_2,q_3...q_N], q_i - 3x1
%
%       p = H_q2p(1:3,:)*[q(:,idx),ones(1,N)]
%
%   [idx,err] = FINDPOINTCORRESPONDENCE(q,p) also returns the summed least
%   squares error (see pointsToSE3).
%
%   Note: NaN values in any point will assume that that point was occluded 
%         during measurements and it will be removed from the set used to
%         calculte the rigid body motion.
%
%   See also pointsToSE3, invSE
%
%   (c) M. Kutzer 03Jan2017, USNA

%% Check Inputs
narginchk(2,2);

N = size(p,2);
if size(q,2) ~= N
    error('p and q must be the same dimension.');
end

if size(q,1) ~= 3
    error('The first input must be a 3xN numeric array.');
end

if size(p,1) ~= 3
    error('The second input must be a 3xN numeric array.');
end

%% Check for occluded points
[~,j] = find(isnan(q));
q(:,j) = [];
p(:,j) = [];
[~,j] = find(isnan(p));
q(:,j) = [];
p(:,j) = [];

N = size(p,2);
if size(p,2) < 4
    warning('There must be at least 4 unoccluded corresponding points to calculated the relative rigid body motion.');
    idx = [];
    err = inf;
    return
end

%% Define all possible index combinations
idx_all = perms(1:N);
M = size(idx_all,1);

idx = cell(1,M);
err = inf(1,M);
for i = 1:M
    idx{i} = idx_all(i,:);
    [~,err(i)] = pointsToSE3(q(:,idx{i}),p);
end
[~,idx_srt] = sort(err,2,'ascend');

idx = idx{ idx_srt(1) };
err = err( idx_srt(1) );
