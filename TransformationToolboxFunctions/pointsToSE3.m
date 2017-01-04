function [H_q2p,err] = pointsToSE3(q,p)
% POINTSTOSE3 finds the best fit rigid body between two sets of point
% correspondences.
%
%   H_q2p = POINTSTOSE3(q,p) This function finds a rigid body  
%   transformation that best relates, in a least squares sense, Frame Q 
%   (the frame to which all points in "q" are referenced) to Frame P (the 
%   frame to which all points in "p" are referenced). The resultant 
%   transformation can be used as follows: 
%
%   Given two sets of N corresponding points:
%       p = [p_1,p_2,p_3...p_N], p_i - 3x1 
%       q = [q_1,q_2,q_3...q_N], q_i - 3x1
%
%   Approximations can be calculated:
%       p_est = H_q2p(1:3,:)*[q; ones(1,N)]
%       
%       H_p2q = invSE(H_q2p);
%       q_est = H_p2q(1:3,:)*[p; ones(1,N)]
%
%   This calculation assumes that there is a correspondence between each 
%   point contained in q and p (i.e. q(:,i) <--> p(:,i)).
%
%   [H_q2p,err] = POINTSTOSE3(q,p) also returns the summed least squares 
%   error:
%   
%   p_est = H_q2p(1:3,:)*[q; ones(1,N)]
%
%   err = sum( sqrt(sum((p - p_est{i}).^2, 1)) );
%
%   Note: NaN values in any point will assume that that point was occluded 
%         during measurements and it will be removed from the set used to
%         calculte the rigid body motion.
%
%   Note: Each of these methods will produce a reflection if the result is
%         better than that of a valid rigid body motion. If the best fit is
%         a reflection for both methods, an empty set is returned with an
%         associated error of "Inf."
%   
%   METHOD 1
%   Special Considerations:
%       1) We are assuming the following about covariance of the error 
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
%   See also findPointCorrespondence, invSE
%
%   (c) M. Kutzer 10July2015, USNA

% Updates
%   03Jan2017 - Updated to include automatic best fit selection.

%TODO - implement special cases 3.2.4 from [2]

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
    H_q2p = [];
    err = inf;
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
    ZERO = 1e-7;    % approximately zero
    d = diag(D);    % get diagonal elements
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

%% Check for valid results
% Display the associated messages from checking for valid rigid body
% transformation
%{
for i = 1:2
    [~,msg] = isSE(H{i});
    fprintf('Method %d: %s\n',i,msg);
end
%}

% Check for valid result for Method 1
if ~isSE(H{1})
    % Remove invalid rigid body transfomation
    %disp(H{1})
    H{1} = [];
end

% Check for valid result for Method 2
if isSE(H{2})
    % Note: The solution as published returns the inverse of the desired 
    %       output.
    H{2} = invSE(H{2});
else
    % Remove invalid rigid body transfomation
    %disp( H{2} )
    H{2} = [];
end

%% Select best solution
p_est = cell(size(H));
delta_p = nan(size(H));
for i = 1:numel(H)
    if ~isempty(H{i})
        % Estimate for p
        p_est{i} = H{i}(1:3,:)*[q; ones(1,N)];
        % Summed least squares error
        delta_p(i) = sum( sqrt(sum((p - p_est{i}).^2, 1)) );
    else
        delta_p(i) = inf;
    end
end

if sum( isfinite(delta_p) ) > 0
    [~,idx] = sort(delta_p,2,'ascend');

    H_q2p = H{ idx(1) };
    err = delta_p( idx(1) );
else
    % warning('No valid solution found. Closest solution may be a reflection.');
    H_q2p = [];
    err = inf;
end
