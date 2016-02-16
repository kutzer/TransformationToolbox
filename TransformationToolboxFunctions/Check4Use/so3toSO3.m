function R = so3toSO3(M)
% so3toSO3 calculates the exponential of a skew-symmetric 3×3 matrix is 
% computed by means of the Rodrigues formula.

%% check inputs
[m,n] = size(M);
if m == 3 && n == 3
    v = vee(M,'fast');
elseif m == 3 && n == 1
    v = M;
    M = wedge(v);
elseif m == 1 && n == 3
    v = transpose(M);
    M = wedge(v);
else
    error('Input must be a 3x3 skew-symmetric matrix or a 3-element vector.')
end

%% calculate Rotation matrix
theta = norm(v);
S = M./theta;

R = eye(3) + sin(theta)*S + (1-cos(theta))*S^2;