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
%                  specified or ZERO = 1e-8, a default value is used.
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
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast

%% Check Inputs
narginchk(1,4);

% Check H
if ~iscell(H)
    error('Input must be defined as an N-element cell array.');
end

% Set defaults
throwErr = false;
ZERO = 1e-8;
METHOD = 'Coupled';

% Parse ZERO and "fast" values
[ZERO,throwErr,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,throwErr);

% Parse variable inputs
for i = 1:numel(cellOut)
    switch lower( class(cellOut{i} ))
        case 'char'
            METHOD = cellOut{i};
        case 'string'
            METHOD = char( cellOut{i} );
        otherwise
            % TODO - check unused elements of cellOut
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

%% Check for valid list of values
[H,info] = validCorrespondenceSE(H,ZERO);

idxRmv = find(info.RemoveBin(1,:));
idxAlt = find(info.AlteredBin(1,:));
msg = '';
if ~isempty(idxRmv)
    msg = sprintf('%sThe following values of H were removed because they are invalid:\n\t[ ',msg);
    idxStr = sprintf('%d ',idxRmv);
    msg = sprintf('%s%s]\n\tSee validCorrespondenceSE for more info.',msg,idxStr);
end

if ~isempty(idxAlt)
    if ~isempty(msg)
        msg = sprintf('%s\n\n',msg);
    end
    msg = sprintf('%sThe following values of H were altered using nearestSE:\n\t[ ',msg);
    idxStr = sprintf('%d ',idxAlt);
    msg = sprintf('%s%s]\n\tSee validCorrespondenceSE for more info.',msg,idxStr);
end

if ~isempty(msg)
    if throwErr
        error(msg);
    else
        warning(msg);
    end
end

%% Number of elements of H
N = numel(H);           

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