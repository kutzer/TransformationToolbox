function muH = meanSE(H,fast)
% MEANSE approximates the mean of a sample of rigid body transformations.
%   muH = MEANSE(H) approximates the mean of a sample of rigid body
%   transformations defined as elements of an N-element cell array, H.
%
%   Input(s)
%       H - n-element cell array containing rigid body transformations
%           H{i} \in SE(M) - ith sample
%
%   Output(s)
%       muH - (M+1)x(M+1) array containing the mean of transformations
%           muH \in SE(M)
%       H{i} \in SE(M) - ith sample
%
%   References:
%   [1] A.W. Long, K.C. Wolfe, M.J. Mashner, & G.S. Chirikjian, "The Banana
%       Distribution is Gaussian: A Localization Study with Exponential
%       Coordinates." Robotics: Science and Systems VIII (2013): 265.
%
%   See also covSE, isSE
%
%   M. Kutzer, 04Jan2017, USNA

% Updates
%   19Apr2021 - Replaced inv with invSE
%   20Apr2021 - Updated to include "fast"
%   02Nov2021 - Updated to use "fast" invSE
%   02Nov2021 - Updated to force real skew-symmetry in delta_h

%% Check Inputs
narginchk(1,2);
if nargin < 2
    fast = false;
end

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
    msg = sprintf('One or more elements of your sample are not valid elements of SE:\n');
    for i = 1:numel(msgs)
        msg = [msg, sprintf('\tElement %d: %s\n',idx(i),msgs{i})];
    end
    if ~fast
        error(msg);
    else
        fprintf('%s\n',msg);
    end
end

%% Calculate mean
N = numel(H);   % Total number of samples
M = size(H,1);  % Dimension of matrix

% Initialize mean
muH{1} = eye(M);    % Current mean
muH{2} = inf(M);    % Updated mean

goFlag = true;
while goFlag
    summation = zeros(M);
    for i = 1:N
        delta_h = logm( invSE(muH{1},1) * H{i} );
        
        % Force skew-symmetry in rotation portion of screw
        %   delta_h \in se(3)
        %   delta_r \in so(3)
        %   ~ Forcing delta_r to be skew-symmetric should account for 
        %     round-off errors
        %   TODO - find documentation reinforcing this approach
        delta_r = delta_h(1:3,1:3);
        delta_r = forceRealSkewSymmetric(delta_r);
        delta_h(1:3,1:3) = delta_r;
        
        summation = summation + delta_h;
    end
    muH{2} = muH{1} * expm( (1/N)*summation );
    
    if isZero(muH{2} - muH{1})
        goFlag = false;
        muH = muH{2};
        break
    else
        muH{1} = muH{2};
    end
end