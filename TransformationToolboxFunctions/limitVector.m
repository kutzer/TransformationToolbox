function [v_out,ratio,inLimits] = limitVector(v,v_lims,method)
% LIMITVECTOR scales the magnitude of a vector to force the infinity norm
% to a value that does not exceed specified limits.
%   v_out = limitVector(v,v_lims)
%
%   Input(s)
%       v      - Nx1 array defining an N-dimensional vector
%       v_lims - Nx2 array specifying the lower and upper limits for each
%                dimension of the vector vector
%           vlims(1,:) = [dimension 1 lower limit, dimension 1 upper limit]
%           vlims(2,:) = [dimension 2 lower limit, dimension 2 upper limit]
%           ...        = ...
%           vlims(N,:) = [dimension N lower limit, dimension N upper limit]
%
%           NOTE: Lower limits must be negative, and upper limits must be
%                 positive in the current version of this code.
%
%       method - [OPTIONAL] character array specifying method
%           ['Inf']       - Use infinity norm and preserve direction (this
%                           is the default method)
%           'ElementWise' - Use an element-wise limit (direction is not
%                           preserved)
%
%   Output(s)
%       v_out - Nx1 array defining an N-dimensional vector bounded by the
%               specified limits such that:

%           'Inf' method
%               (1) Direction is preserved
%                   v./norm(v) = v_out./norm(v_out)
%               (2) New vector is within specified limits
%                   v_out >= v_lims(:,1) & v_out <= v_lims(:,2)
%           'ElementWise' method
%               (1) All elements outside of a specified limit are set to
%                   that limit
%
%       ratio - scalar constant relating v_out to v (applies only to 'Inf'
%               method, otherwise ratio = [])
%
%               v_out = v./ratio
%
%   M. Kutzer, 28Feb2022, USNA

debug = false;

%% Check input(s)
narginchk(2,3);

if nargin < 3
    method = 'Inf';
end

N = numel(v);
v = reshape(v,N,1);
if size(v_lims,1) ~= N || size(v_lims,2) ~= 2
    error('Specified lower/upper limits must be %dx2.',N);
end

bin = v_lims(:,1) < 0;
if nnz(bin) ~= N
    error('Lower limits must be specified as less than 0.')
end

bin = v_lims(:,2) > 0;
if nnz(bin) ~= N
    error('Lower limits must be specified as less than 0.')
end

%% Limit the vector
% Check if vector is within the designated limits
ZERO = 1e-6;
if nnz( (v - v_lims(:,1)) < -ZERO ) == 0 && nnz( (v - v_lims(:,2)) >  ZERO ) == 0
    % Vector is within specified limits
    v_out = v;
    switch lower( method )
        case 'inf'
            ratio = 1;
        otherwise
            ratio = [];
    end
    inLimits = true;
    return
end

% Vector is outside specified limits
inLimits = false;
% Bound desired vector using limits
switch lower( method )
    case 'inf'
        % Define maximum ratio
        ratios = [v,v]./v_lims;
        ratio = max(max( abs(ratios) ));

        % Preserve vectors within specified the limits
        %   (For ratio < 1, ratio = 1)
        if ratio < 1
            ratio = 1;
        end

        % Define vector within specified limits
        v_out = v./ratio;
    case 'elementwise'
        tf_min = v < v_lims(:,1);
        tf_max = v > v_lims(:,2);

        % Define vector with specified limits
        v_out = v;
        v_out(tf_min) = v_lims(tf_min,1);
        v_out(tf_max) = v_lims(tf_max,2);

        % Define ratio
        ratio = [];
    otherwise
        error('Method "%s" is not valid.',method);
end

%% Check properties (debugging)
if debug
    fprintf('(1) Direction is preserved\n');
    fprintf('\t\\hat{v    } = [');
    fprintf('%f ',v./norm(v));
    fprintf(']\n');
    fprintf('\t\\hat{v_out} = [');
    fprintf('%f ',v_out./norm(v_out));
    fprintf(']\n');
    fprintf('\n');

    fprintf('(2) New vector is within specified limits\n');
    fprintf('\tv_out >= v_lims(:,1) = [');
    fprintf('%d ',v_out >= v_lims(:,1));
    fprintf(']\n');
    fprintf('\tv_out <= v_lims(:,2) = [');
    fprintf('%d ',v_out <= v_lims(:,2));
    fprintf(']\n');

    fprintf('(*) Output vector and bounds\n');
    fprintf('\tv_lims(:,1) = [');
    fprintf('%f ',v_lims(:,1));
    fprintf(']\n');
    fprintf('\t      v_out = [');
    fprintf('%f ',v_out);
    fprintf(']\n');
    fprintf('\tv_lims(:,2) = [');
    fprintf('%f ',v_lims(:,2));
    fprintf(']\n');
    fprintf('\n');
end