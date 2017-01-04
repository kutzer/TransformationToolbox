function [SigmaH,muH] = covSE(varargin)
% COVSE approximates the covariance of a sample of rigid body
% transformations.
%   SigmaH = COVSE(H) approximates the covariance of a sample of rigid body
%   transformations defined as elements of an N-element cell array, H.
%
%       H{i} \in SE(M) - ith sample
%   
%   SigmaH = COVSE(H,muH) approximates the covariance of a sample of rigid
%   body transformations defined as elements of an N-element cell array, H
%   given the mean muH.
%
%   [SigmaH,muH] = COVSE(...) returns both the covariance matrix and the
%   mean (either provided or calculated depending on the provided inputs).
%
%   References:
%   [1] A.W. Long, K.C. Wolfe, M.J. Mashner, & G.S. Chirikjian, "The Banana
%       Distribution is Gaussian: A Localization Study with Exponential 
%       Coordinates." Robotics: Science and Systems VIII (2013): 265.
%
%   See also meanSE, isSE
%
%   M. Kutzer, 04Jan2017, USNA

%% Check Inputs
narginchk(1,2);

H = varargin{1};

if ~iscell(H)
    error('Input must be defined as an N-element cell array.');
end

msgs = {};
idx = [];
N = numel(H);
for i = 1:N
    [bin,msg] = isSE(H{i});
    if ~bin
        msgs{end+1} = msg;
        idx(end+1) = i;
    end
end

if ~isempty(idx)
    errmsg = sprintf('One or more elements of your sample are not valid elements of SE:\n');
    for i = 1:numel(msgs)
        errmsg = [errmsg, sprintf('\tElement %d: %s\n',idx(i),msgs{i})];
    end
    error(errmsg);
end

if nargin > 1
    % Check provided muH
    muH = varargin{2};
    [bin,msg] = isSE(muH);
    if ~bin
        errmsg = sprintf('Provided mean is not a valid element of SE:\n\t%s\n',msg);
        error(errmsg);
    end
    
    if size(muH,1) ~= size(H{1},1)
        error('Provided mean must be the same dimension as the provided samples.');
    end
else
    % Calculate muH
    muH = meanSE(H);
end
        
%% Calculate covariance
N = numel(H);
M = numel( veeSE(logm(muH),'fast') );
SigmaH = zeros(M);
for i = 1:N
    y_i = veeSE(logm( inv(muH)*H{i} ),'fast');
    SigmaH = SigmaH + y_i*transpose(y_i);
end
SigmaH = (1/N) * SigmaH;