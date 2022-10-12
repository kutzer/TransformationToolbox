function [v_close,v_sort,nv_sort,idx_sort] = findClosestVector(v_all,v)
% FINDCLOSESTVECTOR finds the vector in a set closest to a reference 
% vector. Closest is defined in terms of the Euclidean norm.
%   v_close = FINDCLOSESTVECTOR(v_all,v)
%
%   [v_close, v_sort] = FINDCLOSESTVECTOR(___)
%
%   [v_close,v_sort,nv_sort] = FINDCLOSESTVECTOR(___)
%
%   [v_close,v_sort,nv_sort,idx_sort] = FINDCLOSESTVECTOR(___)
%
%   Input(s)
%       v_all   - NxM array containing M Nx1 vectors
%       v       - Nx1 vector
%
%   Output(s)
%       v_close - Nx1 vector closest to v (per the Euclidean norm). If 
%                 multiple vectors are found with the same overall distance
%                 from the desired vector, the first is returned.
%       v_sort  - NxM array containing M Nx1 vectors sorted by Euclidean
%                 distance (shortest to closest)
%       nv_sort - 1xM array providing the Euclidean norm of the difference
%                 |v_all - v|
%       idx_sort - 1xM array containing the sorting index such that
%                  v_sort = v_all(:,idx_sort)
%
%   M. Kutzer, 15Dec2016, USNA

% Updates:
%   12Oct2022 - Added idx_sort as optional output and updated documentation

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
[nv_sort,idx] = sort(nv);

v_close = v_all(:,idx(1));
v_sort  = v_all(:,idx);
idx_sort = idx;
