function [SigmaH,muH] = covSE(H,varargin)
% COVSE approximates the covariance of a sample of rigid body
% transformations.
%   SigmaH = COVSE(H) approximates the covariance of a sample of rigid body
%   transformations defined as elements of an N-element cell array, H.
%
%   SigmaH = COVSE(___,muH) allows the user to specify the mean
%
%   SigmaH = MEANSE(__,throwErr)
%
%   SigmaH = MEANSE(__,ZERO)
%
%   SigmaH = MEANSE(__,METHOD)
%
%   [SigmaH,muH] = COVSE(___) returns both the covariance matrix and the
%   mean (either provided or calculated depending on the provided inputs).
%
%   Input(s)
%         H      - n-element cell array containing rigid body transformations
%                  H{i} \in SE(M) - ith sample
%       muH      - [OPTIONAL] (M+1)x(M+1) array containing the mean of 
%                  transformations
%                  muH \in SE(M)
%       throwErr - [OPTIONAL] true/false logical value indicating whether 
%                  or not an error should be thrown if one or more elements 
%                  of H are not valid elements of SE(M).
%                    throwErr = true    - Throw an error if H{i} 
%                                         \notin SE(M)
%                    throwErr = [false] - Warn, but do not throw an error
%       ZERO     - [OPTIONAL] positive value that is sufficiently close to 
%                  zero or assumed zero (e.g. ZERO = 1e-8). If ZERO is not  
%                  specified or ZERO = [], a default value is used.
%       METHOD   - [OPTIONAL] character array specifying method for 
%                  calculating the covariance.
%                     'Coupled' - [DEFAULT] follows [1] finding the mean of 
%                                 the banana distribution
%                   'Decoupled' - finds the rotational covariance 
%                                 independent of the position covariance 
%
%   Output(s)
%       SigmaH - 
%          muH - (M+1)x(M+1) array containing the mean of transformations
%                 muH \in SE(M)
%
%   References:
%   [1] A.W. Long, K.C. Wolfe, M.J. Mashner, & G.S. Chirikjian, "The Banana
%       Distribution is Gaussian: A Localization Study with Exponential 
%       Coordinates." Robotics: Science and Systems VIII (2013): 265.
%
%   See also meanSE, isSE
%
%   M. Kutzer, 04Jan2017, USNA

% TODO - match the syntax of meanSE

%% Check Inputs
narginchk(1,5);

% Check H
if ~iscell(H)
    error('Input must be defined as an N-element cell array.');
end

% Set defaults
throwErr = false;
ZERO = [];
METHOD = 'Coupled';
muH = [];

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
                if ismatrix(varargin{i}) && numel(varargin{i}) == numel(H{1})
                    muH = varargin{i};
                elseif numel(varargin{i}) == 1
                    if varargin{i} == 0 || varargin{i} == 1
                        throwErr = logical(varargin{i});
                    else
                        ZERO = varargin{i};
                    end
                else
                    error('Numeric optional inputs must be scalar values or a valid element of SE(M).');
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

%% Calculate mean (if not provided)
if isempty(muH)
    muH = meanSE(H,throwErr,ZERO,METHOD);
end

% TODO - check mean

%% Check number of elements of H
if N < 0
    error('No valid elements of H present. Unable to calculate mean.');
end

% Initialize covariance
M = numel( seBasis( dims(1)-1 ) );
SigmaH = zeros(M);

if N == 1
    fprintf('\nOne valid element of H, returning zero covariance value.\n');
    return
end


%% Calculate covariance
fast = true;
for i = 1:N
    try
        y_i = veeSE(logSE(invSE(muH,fast)*H{i},fast),fast);
    catch 
        invSE(muH)*H{i}
        break
    end
    SigmaH = SigmaH + y_i*transpose(y_i);
end
SigmaH = (1/N) * SigmaH;