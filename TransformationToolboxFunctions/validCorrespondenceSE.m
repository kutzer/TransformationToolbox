function varargout = validCorrespondenceSE(varargin)
% VALIDCORRESPONDENCESE checks all elements of cell arrays containing
% elements of SE(N) for valid properties of SE(N), and corresponding values
% of SE(N). Bad values are removed and valid corresponence arrays are
% returned.
%   [H_a2b_v,H_c2d_v,...,info] = VALIDCORRESPONDENCESE(H_a2b,H_c2d,...)
%
%   [H_a2b_v,H_c2d_v,...,info] = VALIDCORRESPONDENCESE(___,ZERO)
%
%   Input(s)
%       H_a2b - k-element cell array containing approximate values of SE(N)
%       H_c2d - k-element cell array containing approximate values of SE(N)
%       ...
%       H_y2z - k-element cell array containing approximate values of SE(N)
%       ZERO - positive value that is sufficiently close to zero to be
%              assumed zero (e.g. ZERO = 1e-8). If ZERO is not specified or
%              if ZERO = [], a default value is used.
%
%   Output(s)
%       H_a2b_v - m-element cell array containing values of SE(N)
%       H_c2d_v - m-element cell array containing values of SE(N)
%       ...
%       H_y2z_v - m-element cell array containing values of SE(N) 
%       info    - structured array describing removed/altered transforms
%           info.RemoveMsg - n x k cell array containing descriptions of
%                            why transformations were removed 
%           info.RemoveBin - n x k logical array defining removed
%                            transformations
%           info.AlteredBin - n x k logical array defining altered
%                             transformations
%
%           NOTE: n is the total number of cell arrays provided
%
%   M. Kutzer, 08Sep2022, USNA

% Update(s)
%   09Sep2022 - Updated documentation, added info.AlteredBin, and corrected
%               errors

%% Set default values
ZERO = [];

%% Parse Inputs
n = numel(varargin); % Note we can use "nargin" here instead
tf_H = false(1,n);

for i = 1:numel(varargin)
    if iscell(varargin{i})
        tf_H(i) = true;
    end
end

if nnz(~tf_H) > 1
    error('Only one non-cell input can be specified.')
end

if nnz(~tf_H) > 0
    ZERO = varargin{~tf_H}; % Isolate zero value
    varargin(~tf_H) = [];   % Isolate candidate transforms
end
idx_H = find(tf_H);     % Find index values of varargin to use inputname.m

% Check zero
if ~isempty(ZERO)
    if numel(ZERO) > 1
        error('ZERO value must be a scalar value.');
    end

    if ZERO < 0
        error('ZERO value must be greater or equal to 0.')
    end
end

%% Check number of elements of each cell
n = numel(varargin);    % Number of corresponding sets
m = zeros(1,n);         % Number of elements in each set
for i = 1:n
    H_i = varargin{i};
    m(i) = numel(H_i);
end

m_star = mode(m); % <--- most frequently occuring value
tf_Bad = m ~= m_star;
if nnz( tf_Bad ) > 0
    % Initialize message
    msg = sprintf('Cell arrays containing corresponding matrices must be the same length.\n');

    for i = 1:n
        if tf_Bad(i)
            % Get ith input name
            str_i = inputname(idx_H(i));
            % Update message
            msg = sprintf('%s\tnumel(%s) = %d ~= %d\n',msg,str_i,m(i),m_star);
        end
    end
    error(msg);
end

m = m_star; % Update number of elemens in each set

%% Check dimensions of matrix elements
tfRemove  = false(n,m); % Elements to remove
tfAltered = false(n,m); % Elements altered, but not removed
tfMsg = cell(n,m);      % Details on why element was removed
mm = zeros(n,m);        % Array of 1st dimension value

for i = 1:n
    for j = 1:m
        H_ij = varargin{i}{j};  % Should be \in SE(N)
        
        % Check if element is a matrix
        if ~ismatrix(H_ij)
            % Get ith input name
            str_i = inputname(idx_H(i));
            % Define message
            msg = sprintf('%s{%d} is not a matrix',str_i,j);
            tfMsg{i,j} = msg;
            % Update remove list
            tfRemove(i,j) = true;
            continue
        end
        
        % Check if element is square
        [mm(i,j),nn] = size(H_ij);
        if mm(i,j) ~= nn
            % Get ith input name
            str_i = inputname(idx_H(i));
            % Define message
            msg = sprintf('%s{%d} is not a square matrix: %d x%d',str_i,j,...
                mm(i,j),nn);
            tfMsg{i,j} = msg;
            % Update remove list
            tfRemove(i,j) = true;
            continue
        end
        
        % Check if dimension meets minimum requirement
        if nn < 3
            % Get ith input name
            str_i = inputname(idx_H(i));
            % Define message
            msg = sprintf('%s{%d} is %d x %d, must be at least 3 x 3',str_i,j,...
                mm(i,j),nn);
            tfMsg{i,j} = msg;
            % Update remove list
            tfRemove(i,j) = true;
            continue
        end

    end
end

%% Check for common dimensions
mm_star = mode( mm(~tfRemove) );        % Most common matrix dimension
tfBadDim = (mm ~= mm_star) & ~tfRemove; % New bad values
[i_all,j_all] = find(tfBadDim);
for k = 1:numel(i_all)
    i = i_all(k);
    j = j_all(k);
    % Get ith input name
    str_i = inputname(idx_H(i));
    % Define message
    msg = sprintf('%s{%d} is %d x %d, should be %d x %d',str_i,j,...
        mm(i,j),mm(i,j),mm_star,mm_star);
    tfMsg{i,j} = msg;
    % Update remove list
    tfRemove(i,j) = true;
end

%% Remove all correspondences that have at least one bad value
tfBadCor = repmat( sum(tfRemove,1) > 0, n, 1 );
[i_all,j_all] = find(tfBadCor);
for k = 1:numel(i_all)
    i = i_all(k);
    j = j_all(k);
    if tfBadCor(i,j) && ~tfRemove(i,j)
        % Get ith input name
        str_i = inputname(idx_H(i));
        % Define message
        msg = sprintf('%s{%d} corresponds to a bad value',str_i,j);
        tfMsg{i,j} = msg;
    end
end

tfRemove = tfBadCor;

%% Map to nearest
II = eye(mm_star);               % Identifty matrix
[i_all,j_all] = find(~tfRemove);
for k = 1:numel(i_all)
    i = i_all(k);
    j = j_all(k);

    H_ij = varargin{i}{j};      % Should be \in SE(N)
    H_ij_SE = nearestSE(H_ij);  % Force valid \in SE(N)

    if ~isSE(H_ij_SE,ZERO)
        % Get ith input name
        str_i = inputname(idx_H(i));
        % Define message
        msg = sprintf('%s{%d} does not have a valid nearest element of SE(%d)',str_i,j,mm_star-1);
        tfMsg{i,j} = msg;
        % Update remove list
        tfRemove(i,j) = true;
        continue
    end

    if ~isZero(H_ij*invSE(H_ij_SE,true)-II, ZERO)
        % Get ith input name
        str_i = inputname(idx_H(i));
        % Define message
        msg = sprintf('%s{%d} has a larger than ZERO change when mapped to nearest SE(%d)',str_i,j,mm_star-1);
        tfMsg{i,j} = msg;
        tfAltered(i,j) = true;

        % ??DO NOT REMOVE??
    end

    % Update with nearest value
    varargin{i}{j} = H_ij_SE;
end

%% Remove all correspondences that have at least one bad value
tfBadCor = repmat( sum(tfRemove,1) > 0, n, 1 );
[i_all,j_all] = find(tfBadCor);
for k = 1:numel(i_all)
    i = i_all(k);
    j = j_all(k);
    if tfBadCor(i,j) && ~tfRemove(i,j)
        % Get ith input name
        str_i = inputname(idx_H(i));
        % Define message
        msg = sprintf('%s{%d} corresponds to a bad value',str_i,j);
        tfMsg{i,j} = msg;
    end
end

tfRemove = tfBadCor;

%% Package outputs
for i = 1:nargout
    if i <= n
        varargout{i} = varargin{i}(~tfRemove(i,:));
    elseif i == n+1
        info.RemoveMsg  = tfMsg;
        info.RemoveBin  = tfRemove;
        info.AlteredBin = tfAltered;
        varargout{i} = info;
    else
        warning('Maximum valid outputs is %d, output %d is empty.',n+1,i);
        varargout{i} = [];
    end
end