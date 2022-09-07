function muH = meanSE(H,varargin)
% MEANSE approximates the mean of a sample of rigid body transformations.
%   muH = MEANSE(H) approximates the mean of a sample of rigid body
%   transformations defined as elements of an N-element cell array, H.
%
%   muH = MEANSE(__,throwErr)
%
%   muH = MEANSE(__,ZERO)
%
%   muH = MEANSE(__,METHOD)
%
%   Input(s)
%       H        - n-element cell array containing (N+1)x(N+1) array
%                  elements of SE(N)
%                  H{i} \in SE(N) - ith sample
%       throwErr - [OPTIONAL] true/false logical value indicating whether 
%                  or not an error should be thrown if one or more elements 
%                  of H are not valid elements of SE(N).
%                    throwErr = true    - Throw an error if H{i} 
%                                         \notin SE(N)
%                    throwErr = [false] - Warn, but do not throw an error
%       ZERO     - [OPTIONAL] positive value that is sufficiently close to 
%                  zero or assumed zero (e.g. ZERO = 1e-8). If ZERO is not  
%                  specified or ZERO = [], a default value is used.
%       METHOD   - [OPTIONAL] character array specifying method for 
%                  calculating the mean.
%                     'Coupled' - [DEFAULT] follows [1] finding the mean of 
%                                 the banana distribution
%                   'Decoupled' - finds the rotational mean independent of
%                                 the position mean 
%
%   Output(s)
%       muH - (N+1)x(N+1) array containing the mean of transformations such
%             that muH \in SE(N)
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
%   18Nov2021 - Replaced "fast" with "throwErr"
%   18Nov2021 - Updated to include ZERO definition
%   22Nov2021 - Fixed varargin{1} error
%   30Nov2021 - Special case for 2-transformation mean
%   10Mar2022 - Added else for islogical(varargin{1}) 
%   01Sep2022 - Added coupled/decoupled option and revised parsing of
%               optional inputs
%   02Sep2022 - Corrected/removed SE(3) dimension assumptions

%% Check Inputs
narginchk(1,4);

% Check H
if ~iscell(H)
    error('Input must be defined as an N-element cell array.');
end

% Set defaults
throwErr = false;
ZERO = [];
METHOD = 'Coupled';

% Parse variable inputs
if nargin > 1
    for i = 1:numel(varargin)
        switch lower( class(varargin{i} ))
            case 'char'
                METHOD = varargin{i};
            case 'string'
                METHOD = char( varargin{i} );
            case 'logical'
                throwErr = varargin{i};
            otherwise
                if numel(varargin{i}) ~= 1
                    error('Numeric optional inputs must be scalar values.');
                end

                if varargin{i} == 0 || varargin{i} == 1
                    throwErr = logical(varargin{i});
                else
                    ZERO = varargin{i};
                end
        end
    end
end

% Check method
switch lower(METHOD)
    case 'coupled'
        % Good
    case 'decoupled'
        % Good
    otherwise
        error('"%s" is not a valid method. Please use "Coupled" or "Decoupled".',METHOD);
end

% Check zero
if ZERO < 0
    error('ZERO value must be greater or equal to 0.')
end

%% Check values H for valid SE(N)
% Initialize variables
N = numel(H);           % Number of elements of H
dims = zeros(N,1);      % Dimensions of each element of H
msgs = {};              % Messages desribing bad values of H
idx = [];               % Indices associated with bad valus of H
H_inSE = H;             % Corrected values of H
nearSE = true(size(H)); % Logical array describing corrected values  
                        %   resulting in a valid elements of SE(M)
nearSEmsg = cell(N,1);  % Cell array describing bad nearest element of SE(M)

% Check values of H
for i = 1:N
    % Check dimensions of H
    if ~ismatrix(H{i})
        error('All elements of H must be square matrices.');
    end
    dims(i) = size(H{i},1);
    % Check ith element of H
    [bin,msg] = isSE(H{i},ZERO);
    if ~bin
        % Compile messages and add index for error/warning reporting
        msgs{end+1} = msg;
        idx(end+1) = i;

        % Attempt to map bad elements to SE(M)
        if ~throwErr
            H_inSE{i} = nearestSE(H{i});
            [nearSE(i),nearSEmsg{i}] = isSE(H_inSE{i},ZERO);
        end
    end
end

% Throw error for inconsistent or incorrect dimensions
if nnz(dims < 3) > 0
    error('All elements of H must be 3x3 or larger.');
end
if nnz(dims(1) ~= dims) > 0
    error('All elements of H must be the same size.');
end

% Throw error/warning for bad elements of SE(M)
if ~isempty(idx)
    msg = sprintf('One or more elements of your sample are not valid elements of SE:\n');
    for i = 1:numel(msgs)
        msg = [msg, sprintf('\tElement %d: %s\n',idx(i),msgs{i})];
    end
    if throwErr
        error(msg);
    else
        fprintf('%s\n',msg);
    end
end

%% Replace/remove bad nearest elements of SE(M)
if numel(msgs) > 0
    badVals = numel(msgs);       % Number of bad values of H
    remVals = nnz(~nearSE);      % Number of values to remove from H
    newVals = badVals - remVals; % Number of nearestSE replaced values of H
    
    % Build message
    msg = [];
    msg = sprintf('%s\t%d elements of H are not valid elements of SE(M)\n',msg,badVals);
    msg = sprintf('%s\t%d elements of H replaced with nearest SE(M)\n'    ,msg,newVals);
    msg = sprintf('%s\t%d elements of H removed\n'                        ,msg,remVals);
    fprintf('%s',msg);
    
    % Replace H with nearest SE(M) replacemets
    H = H_inSE;
    % Remove values without a valid SE(M) element
    H(~nearSE) = [];

    % Update the number of elements of H
    N = numel(H);
end

%% Check number of elements of H
if N < 0
    error('No valid elements of H present. Unable to calculate mean.');
end
if N == 1
    fprintf('\nOne valid element of H, returning single value.\n');
    muH = H{1};
    return
end

%% Get dimensions of H
% NOTE: From this point on, we will refer to elements of H as SE(M-1)
M = size(H{1},1);  % Dimension of matrix

%% Account for coupled/decoupled method
switch lower(METHOD)
    case 'coupled'
        % Do nothing
    case 'decoupled'
        % Isolate rotation and translation
        D = zeros(M-1,N);
        for i = 1:N
            D(:,i) = H{i}(1:(M-1),M);   % Isolate translation
            H{i}(1:(M-1),M) = 0;        % Isolate rotation
        end
        % Calculate translation mean
        muD = mean(D,2);
end

%% Special case (2 transforms)
if numel(H) == 2
    H_1to2 = invSE(H{2})*H{1};
    delta_h = logSE( H_1to2 );
    
    % Force skew-symmetry in rotation portion of screw
    %   delta_h \in se(M-1)
    %   delta_r \in so(M-1)
    %   ~ Forcing delta_r to be skew-symmetric should account for
    %     round-off errors
    %   TODO - find documentation reinforcing this approach
    delta_r = delta_h(1:(M-1),1:(M-1));
    delta_r = forceRealSkewSymmetric(delta_r);
    delta_h(1:(M-1),1:(M-1)) = delta_r;
        
    muH = H{2} * expm( (1/2)*delta_h );
    
    % Apply decoupled method
    switch lower(METHOD)
        case 'decoupled'
            muH(1:(M-1),M) = muD;
    end

    return
end

% %% Special case (2 transforms)
% if numel(H) == 2
%     H_1to2 = invSE(H{2})*H{1};
%     X_1to2 = H_1to2(1:(M-1),M);
%     R_1to2 = H_1to2(1:(M-1),1:(M-1));
%     r_1to2 = vee(logm(R_1to2));
%     dR = expm( wedge( r_1to2./2 ) );
%     dX = X_1to2./2;
%     dH = eye(M);
%     dH(1:(M-1),1:(M-1)) = dR;
%     dH(1:(M-1),M) = dX;
%     muH = H{2} * dH;
%     return
% end

%% Calculate mean

% Initialize mean
muH{1} = eye(M);    % Current mean
muH{2} = inf(M);    % Updated mean

goFlag = true;
while goFlag
    summation = zeros(M);
    for i = 1:N
        delta_h = logSE( invSE(muH{1},1) * H{i} );
        
        % Force skew-symmetry in rotation portion of screw
        %   delta_h \in se(M-1)
        %   delta_r \in so(M-1)
        %   ~ Forcing delta_r to be skew-symmetric should account for 
        %     round-off errors
        %   TODO - find documentation reinforcing this approach
        delta_r = delta_h(1:(M-1),1:(M-1));
        delta_r = forceRealSkewSymmetric(delta_r);
        delta_h(1:(M-1),1:(M-1)) = delta_r;
        
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

% Apply decoupled method
switch lower(METHOD)
    case 'decoupled'
        muH(1:(M-1),M) = muD;
end