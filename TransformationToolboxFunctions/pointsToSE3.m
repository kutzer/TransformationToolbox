function H = pointsToSE3(q,p)
% POINTSTOSE3 finds the best fit rigid body between two sets of point
% correspondences.
%
%   POINTSTOSE3(q,p) This function finds a rigid body motion that best 
%   moves p to q assuming a correspondence between each vector contained 
%   in p and q (i.e. q(:,i) <--> p(:,i)).
%
%   p = H(1:3,:)*[q;ones(1,N)]
%
%   p = [p_1,p_2,p_3...p_N], p_i - 3x1 
%   q = [q_1,q_2,q_3...q_N], q_i - 3x1
% 
%   METHOD 1
%   Special Considerations:
%       1) NaN values in any point will assume that that point was occluded 
%       during measurements and it will be removed from the set used to
%       calculte the rigid body motion.
%
%       2) We are assuming the following about covariance of the error 
%       vectors:
%           W = Wi for all i
%           W = w*eye(3)
%
%   References
%       [1] G.S. Chirikjian & A.B. Kyatkin, "Engineering Applications of 
%       Noncommutative Harmonic Analysis", pp. 472-473, 2001.
%
%   METHOD 2
%   References
%       [2] D.W. Eggert1, A. Lorusso2, R.B. Fisher3, "Estimating 3-D rigid 
%       body transformations: a comparison of four major algorithms," 
%       Machine Vision and Applications (1997) 9: 272–290
%
%   See also 
%
%   (c) M. Kutzer 10July2015, USNA

%TODO - implement special cases 3.2.4 from [2]
%TODO - select best method

%% Check Inputs
narginchk(2,3);

N = size(p,2);
if size(q,2) ~= N
    error('p and q must be the same dimension.');
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
    H = [];
    return
end

%% Calculate relative rigid body motion (METHOD 1)
p_cm = (1/N)*sum(p,2);
q_cm = (1/N)*sum(q,2);

p_rel = bsxfun(@minus,p,p_cm);
q_rel = bsxfun(@minus,q,q_cm);

C = p_rel*transpose(q_rel);

%a_opt = p_rel - ( C*(C'*C)^(-1/2) )*q_rel
a_opt = p_cm - ( C*(transpose(C)*C)^(-1/2) )*q_cm;
A = ( C*(transpose(C)*C)^(-1/2) );

H{1} = eye(4);
H{1}(1:3,1:3) = A;
H{1}(1:3,4) = a_opt;

%% Calculate relative rigid body motion (METHOD 2.1)
p_cm = (1/N)*sum(p,2);
q_cm = (1/N)*sum(q,2);

p_rel = bsxfun(@minus,p,p_cm);
q_rel = bsxfun(@minus,q,q_cm);

C = p_rel*transpose(q_rel);

[U,D,V] = svd(C);

R = V*transpose(U);

if det(R) < 0 % account for reflections
    %TODO - confirm that this step of finding the location of the singular
    %value of C is needed, or if V_prime = [v1,v2,-v3] always
    ZERO = 1e-7; % close enough to zero
    d = undiag(D);
    bin = (abs(d) < ZERO);
    sgn = ones(1,3);
    sgn(bin) = -1;
    V_prime = [sgn(1)*V(:,1),sgn(2)*V(:,2),sgn(3)*V(:,3)];
    R = V_prime*transpose(U);
end

T = q_cm - R*p_cm;

H{2} = eye(4);
H{2}(1:3,1:3) = R;
H{2}(1:3,4) = T;