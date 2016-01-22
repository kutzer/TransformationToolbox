function v = vee(M,options)
%vee converts a 3x3 or 2x2 skew-symmetric matrix into a vector
%   vee(M) converts a 3x3 skew-symmetric matrix "M" into a 3x1 vector "v", 
%   or a 2x2 skew-symmetric matrix "M" into a scalar "v".
%
%   See also wedge isSkewSymmetric
%
%   (c) M. Kutzer 10Oct2014, USNA

%% Default options
%TODO - document "options" 
if nargin < 2
    options = '';
end

%% Check for custom functions
%TODO - check for isSkewSymmetric.m

%% Check M
switch lower(options)
    case 'fast'
        % Do not check for skew symmetry
    otherwise
        if ~isSkewSymmetric(M)
            error('"M" must be skew-symmetric.');
        end
end

%% Calculate v
if size(M,1) == 3
    v(1,1) = M(3,2);
    v(2,1) = M(1,3);
    v(3,1) = M(2,1);
    return
end

if size(M,1) == 2
    v(1,1) = M(2,1);
    return
end

%% $M \in so(N)$
n = size(M,1);
e = soBasis(n);
m = numel(e);
for idx = 1:m
    [i,j] = find(e{idx} == 1);
    v(idx,1) = M(i,j);
end

end