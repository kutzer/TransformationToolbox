function muH = meanSE(H,varargin)
% MEANSE approximates the mean of a sample of rigid body transformations.
%   muH = MEANSE(H) approximates the mean of a sample of rigid body
%   transformations defined as elements of an N-element cell array, H.
%
%   muH = MEANSE(H,tf)
%
%   muH = MEANSE(H,ZERO)
%
%   muH = MEANSE(H,tf,ZERO)
%
%   Input(s)
%       H    - n-element cell array containing rigid body transformations
%              H{i} \in SE(M) - ith sample
%       tf   - true/false value indicating whether or not an error should
%              be thrown if one or more elements of H are not valid 
%              elements of SE(N).
%                   tf = [true] - Throw an error if H{i} \notin SE(N)
%                   tf = false  - Warn, but do not throw an error
%       ZERO - positive value that is sufficiently close to zero to be
%              assumed zero (e.g. ZERO = 1e-8). If ZERO is not specified or
%              if ZERO = [], a default value is used.
%
%   Output(s)
%       muH - (M+1)x(M+1) array containing the mean of transformations such
%             that muH \in SE(M)
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
%   18Nov2021 - Replaced "fast" with "tf"
%   18Nov2021 - Updated to include ZERO definition
%   22Nov2021 - Fixed varargin{1} error
%   30Nov2021 - Special case for 2-transformation mean

%% Check Inputs
narginchk(1,3);

if nargin < 2
    tf = false;
    ZERO = [];
end

if nargin == 2
    % muH = MEANSE(H,tf)
    % muH = MEANSE(H,ZERO)
    if ~islogical(varargin{1})
        % TODO - consider warning and/or enforcing use of logical tf
        %        argument
        if varargin{1} == 0 || varargin{1} == 1
            tf = true;
            ZERO = [];
        else
            tf = false;
            ZERO = varargin{1};
        end
    end
end

if nargin == 3
    % muH = MEANSE(H,tf,ZERO)
    tf = varargin{1};
    ZERO = varargin{2};
end

% TODO - check tf value 

if ~iscell(H)
    error('Input must be defined as an N-element cell array.');
end

% Check values contained in H for valid SE(N)
msgs = {};
idx = [];
N = numel(H);
for i = 1:N
    [bin,msg] = isSE(H{i},ZERO);
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
    if ~tf
        error(msg);
    else
        fprintf('%s\n',msg);
    end
end

%% Special case (2 transforms)
if numel(H) == 2
    H_1to2 = invSE(H{2})*H{1};
    delta_h = logm( H_1to2 );
    
    % Force skew-symmetry in rotation portion of screw
    %   delta_h \in se(3)
    %   delta_r \in so(3)
    %   ~ Forcing delta_r to be skew-symmetric should account for
    %     round-off errors
    %   TODO - find documentation reinforcing this approach
    delta_r = delta_h(1:3,1:3);
    delta_r = forceRealSkewSymmetric(delta_r);
    delta_h(1:3,1:3) = delta_r;
        
    muH = H{2} * expm( (1/2)*delta_h );
    
    return
end

% %% Special case (2 transforms)
% if numel(H) == 2
%     H_1to2 = invSE(H{2})*H{1};
%     X_1to2 = H_1to2(1:3,4);
%     R_1to2 = H_1to2(1:3,1:3);
%     r_1to2 = vee(logm(R_1to2));
%     dR = expm( wedge( r_1to2./2 ) );
%     dX = X_1to2./2;
%     dH = eye(4);
%     dH(1:3,1:3) = dR;
%     dH(1:3,4) = dX;
%     muH = H{2} * dH;
%     return
% end

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
    
    if isZero(muH{2} - muH{1},ZERO)
        goFlag = false;
        muH = muH{2};
        break
    else
        muH{1} = muH{2};
    end
end