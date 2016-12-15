function [v_close,v_sort] = findClosestVector(v_all,v)
% FINDCLOSESTVECTOR finds the vector in a set closest to a reference 
% vector. Closest is defined in terms of the Euclidean norm.
%   v_close = FINDCLOSESTVECTOR(v_all,v) finds the vector in a set (v_all)
%   closest to a reference vector (v). If multiple vectors are found with 
%   the same overall distance from the desired vector, the first is 
%   returned.
%       v_all   - NxM array containing M Nx1 vectors
%       v       - Nx1 vector
%       v_close - Nx1 vector closest to v (per the Euclidean norm)
%
%   [v_close, v_sort] = findClosestVector(v_all,v) returns the closest
%   vector and all vectors from the original set, sorted from closest to
%   furthest. 
%
%   M. Kutzer, 15Dec2016, USNA

%% Check inputs
narginchk(2,2);
v = reshape(v,[],1);

[N,M] = size(v_all);
if numel(v) ~= N
    error('The number of elements in the reference vector must match the first dimension of v_all.');
end

if M < 1 
    v_close = [];
    v_sort  = [];
    return
end

if M == 1
    v_close = v;
    v_sort  = v;
    return
end
    
%% Find closest solution
dv = bsxfun(@minus,v_all,v);    % Find difference
nv = sqrt( sum(dv.^2,1) );      % Find norm of difference
[~,idx] = sort(nv);

v_close = v_all(:,idx(1));
v_sort  = v_all(:,idx);
