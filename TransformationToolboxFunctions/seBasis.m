function e = seBasis(n)
% SEBASIS defines the basis elements of the Lie algebra associated with the
% n-dimensional special Euclidean group, SE(n).
%   SEBASIS(n) This function returns the basis elements of se(n) (the Lie
%   algebra associated with SE(n), nometimes referred to as "little se").
%
%   See also soBasis
%
%   M. Kutzer, 04Jan2017, USNA

%% Calculate basis elements of little so
e_so = soBasis(n);

%% Populate little so elements
N = numel(e_so);
M = size(e_so{1},1);

e = cell(1,N + M);
% Append zeros
for i = 1:N
    e{i} = zeros(M+1);
    e{i}(1:M,1:M) = e_so{i};
end

%% Populate remaining elements
for i = 1:M
    e{N+i} = zeros(M+1);
    e{N+i}(i,M+1) = 1;
end
