function [dX,dq,dR] = configurationDistance(q,qs,DH_hndl)
% CONFIGURATIONDISTANCE calculates the summed distance between assigned
% DH frame origins for given robot configurations.
%   [dX,dq,dR] = CONFIGURATIONDISTANCE(q,qs,DH_hndl)
%
%   Input(s)
%       q       - Nx1 array containing current joint configuration
%       qs      - NxM array containing all candidate joint configurations
%       DH_hndl - DH table function handle such that DH_hndl(q) produces a
%                 DH table for joint configuration q.
%
%   Output(s)
%       dX - average Euclidean distance between all DH frame origins
%            dX(i) - distance between q and qs(:,i) configurations
%       dq - distance between joint configurations
%            dq(i) - |qs(:,i) - q|
%       dR - average rotational distance between all DH frames
%            dR(i) - rotational distance between q and qs(:,i)
%
%   See also distanceSO

%% Check input(s)
N = numel(q);
q = reshape(q,N,1);

if size(qs,1) ~= N
    error('First dimension of q and qs must match.');
end

if ~isa(DH_hndl,'function_handle')
    error('DH_hndl must be a valid function handle.');
end

%% Calculate distances
dh_q = DH_hndl(q);
dh_N = size(dh_q,1);
for i = 1:dh_N
    % Calculate ith frame
    H_i2o_q = DHtableToFkin(dh_q(1:i,:));
    % Isolate rotation
    R_i2o_q{i} = H_i2o_q(1:3,1:3);
    % Isolate translation
    X_i2o_q(:,i) = H_i2o_q(1:3,4);
end


M = size(qs,2);
for j = 1:M
    q_j = qs(:,j);
    dh_qj = DH_hndl(q_j);
    
    dX_j = 0;
    dR_j = 0;
    for i = 1:dh_N
        % Calculate ith frame
        H_i2o_qj = DHtableToFkin(dh_qj(1:i,:));
        % Isolate rotation
        R_i2o_qj = H_i2o_qj(1:3,1:3);
        % Isolate translation
        X_i2o_qj = H_i2o_qj(1:3,4);
        
        % Calculate rotational distance
        dR_j = distanceSO(R_i2o_q{i},R_i2o_qj) + dR_j;
        % Calculate Euclidean distance
        dX_j = norm( X_i2o_q(:,i) - X_i2o_qj ) + dX_j;
    end
    
    dR(j) = dR_j/dh_N;
    dX(j) = dX_j/dh_N;
end

%% Calculate joint configuration distances
deltaQ = qs - repmat(q,1,M);
dq = sqrt( sum(deltaQ.^2,1) );