function e = soBasis(n)
% SOBASIS defines the basis elements of the Lie algebra associated with the
% n-dimensional special orthogonal group, SO(n).
%   SOBASIS(n) This function returns the basis elements of so(n) (the Lie 
%   algebra associated with SO(n), sometimes referred to as "little so").
%
%   References:
%   [1] D. Mortari, "On the Rigid Rotation Concept in n-Dimensional 
%       Spaces?"
%   [2] http://en.wikipedia.org/wiki/Rotation_group_SO(3)
%
%   See also seBasis
%
%   M. Kutzer, 26Feb2015, USNA

%TODO - cleanup and comment

%% 
m = n*(n-1)/2;  % total number of basis elements

EE = reshape((n^2):-1:1,n,n)';
EE = triu(EE) - diag( diag(EE) );

%% 
[~,idx] = sort(reshape(EE,[],1));
idx = idx((end-m+1):end);
for k = 1:m
    [i(k),j(k)] = find(EE == EE(idx(k)));
end

%%
for k = 1:m
    e{k} = zeros(n);
    %V(i,j) = (-1)^(i+j);
    e{k}(i(k),j(k)) = (-1)^(i(k)+j(k));  % upper triangular term
    e{k}(j(k),i(k)) = -(-1)^(i(k)+j(k)); % lower triangualr term (skew-symmetric)
end