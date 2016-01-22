function M = wedge(v)
%wedge converts a 3x1 vector or a scalar into skew-symmetric matrix
%
%   wedge(v) calculates a skew-symmetric matrix from elements of v. This
%   function is defined for $v \in \mathbb{R}^3$ or $v \in \mathbb{R}^1$
%
%   See also vee
%
%   (c) M. Kutzer 09Oct2014, USNA

%Updates
%   27Feb2015 - Updated to include N-dimensional wedge (or hat) operator.

%TODO - check for soBasis.m
m = numel(v);
v = reshape(v,m,[]);

% $M \in so(3)$
if m == 3        
    M = [    0, -v(3),  v(2);...
          v(3),     0, -v(1);...
         -v(2),  v(1),     0];
    return
end

% $M \in so(2)$
if m == 1
    M = [   0, -v(1);...
         v(1),     0];
    return
end


% $M \in so(N)$
%m = numel(v);
%m = n*(n-1)/2;
n = (8*m + 1)^(1/2)/2 + 1/2;
if n ~= round(n)
    error('Invalid dimension for input vector')
end
e = soBasis(n);
M = zeros(size(e{1}));
for i = 1:m
    M = M + v(i)*e{i};
end 

end