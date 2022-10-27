function [H_b2a,H_y2x,varargout] = solveFixedTransforms(H_x2a,H_y2b,varargin)
% SOLVEFIXEDTRANSFORMS solves for the fixed transformations associated
% with corresponding data sets.
%   [H_b2a,H_y2x] = SOLVEFIXEDTRANSFORMS(H_x2a,H_y2b)
%
%   [H_b2a,H_y2x] = SOLVEFIXEDTRANSFORMS(___,ZERO)
%
%   [H_b2a,H_y2x] = SOLVEFIXEDTRANSFORMS(___,fast)
%
%   [H_b2a,H_y2x,stats] = SOLVEFIXEDTRANSFORMS(___)
%
%   NOTE: For large data sets, calculating the "stats" output will require
%         additional computation time.
%
%   Input(s)
%       H_x2a - n-element cell array containing unique elements of SE(N)
%           H_x2a{i} - ith element of H_x2a
%       H_y2b - n-element cell array containing unique elements of SE(N)
%           H_y2b{i} - ith element of H_y2b
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not
%              specified, ZERO = 1e-8 is used.
%       fast - [OPTIONAL] true/false logical value indicating whether to
%              skip checking SE(N). Choosing fast = true ignores specified
%              ZERO.
%                fast = true    - Skip checking if H \in SE(N)
%                fast = [false] - Check if H \in SE(N)
%   Output(s)
%       H_b2a - (N+1)x(N+1) array element of SE(N) defining the fixed
%               transformation relating Frame b to Frame a
%       H_y2x - (N+1)x(N+1) array element of SE(N) defining the fixed
%               transformation relating Frame y to Frame x
%       stats - structured array describing error statistics
%           Statistics related to H_b2a:
%               stats.muH_b2b    - error mean of H_bi2bi
%               stats.muH_a2a    - error mean of H_aj2aj
%               stats.SigmaH_b2b - error covariance of H_bi2bi
%               stats.SigmaH_a2a - error covariance of H_aj2aj
%               stats.H_bi2bj    - k-element cell array containing SE(N)
%               stats.H_ai2aj    - k-element cell array containing SE(N)
%           Statistics related to H_y2x:
%               stats.muH_y2y    - error mean of H_yi2yi
%               stats.muH_x2x    - error mean of H_xj2xj
%               stats.SigmaH_y2y - error covariance of H_yi2yi
%               stats.SigmaH_x2x - error covariance of H_xj2xj
%               stats.H_yi2yj    - k-element cell array containing SE(N)
%               stats.H_xi2xj    - k-element cell array containing SE(N)
%
%       NOTE: stats is calculated as follows
%               LHS_H_bi2aj = H_ai2aj{i}*H_bi2ai;
%               RHS_H_bi2aj = H_bj2aj*H_bi2bj{i};
%               H_bi2bi{i} = invSE(LHS_H_bi2aj) * RHS_H_bi2aj;
%               H_aj2aj{i} = LHS_H_bi2aj * invSE(RHS_H_bi2aj);
%
%               LHS_H_yi2xj = H_xi2xj{i}*H_yi2xi;
%               RHS_H_yi2xj = H_yj2xj*H_yi2yj{i};
%               H_yi2yi{i} = invSE(LHS_H_yi2xj) * RHS_H_yi2xj;
%               H_xj2xj{i} = LHS_H_yi2xj * invSE(RHS_H_yi2xj);
%
%
%   Use Requirements:
%       (1) Correspondence - Elements of H_x2a and H_y2b must correspond
%           (i.e. H_x2a{i} must be captured at the same instant as H_y2b{i})
%       (2) Unique Elements - Elements of H_x2a and H_y2b should be unique
%           (i.e. H_x2a{i} ~= H_x2a{j} \forall i ~= j and
%                 H_x2a{i} ~= H_x2a{j} \forall i ~= j)
%       (3) Fixed Transformations - H_b2a and H_y2x must be fixed
%           transformtaions
%
%   Suppressed Output Option:
%       Syntax Example:
%           [H_b2a,~,stats] = solveFixedTransforms(H_x2a,H_y2b)
%       Description:
%           For cases where only one of two fixed transformations is
%           required, this function can use "detectOutputSuppression.m"
%           To enable this feature, download detectOutputSuppression from
%           the MathWorks File Exchange, and add the function location to
%           your current path.
%       Download Information: detectOutputSuppression, A. Danz
%       https://www.mathworks.com/matlabcentral/fileexchange/79218-detectoutputsuppression
%
%   Example(s):
%       % EXAMPLE 1 - MoCap/Manipulator Co-calibration --------------------
%       %   Given:
%       %       H_e2o - n-element cell array containing elements of SE(3)
%       %               defining the manipulator end-effector (Frame e)
%       %               pose relative to the manipulator base frame
%       %               (Frame o).
%       %       H_t2w - n-element cell array containing elements of SE(3)
%       %               defining a MoCap tool frame (Frame t, fixed
%       %               relative to Frame e) relative to the MoCap world
%       %               frame (Frame w).
%       %       H_b2w - n-element cell array containing elements of SE(3)
%       %               defining a MoCap base frame (Frame b, fixed
%       %               relative to Frame o) relative to the MoCap world
%       %               frame (Frame w).
%       %
%       %   Usage Requirements:
%       %       H_t2e - MoCap tool frame (Frame t) pose relative to
%       %               manipulator end-effector frame (Frame e) *MUST BE
%       %               CONSTANT/FIXED*
%       %       H_b2o - MoCap base frame (Frame b) pose relative to
%       %               manipulator base frame (Frame o) *MUST BE
%       %               CONSTANT/FIXED*
%       %
%       %   NOTE: This example does not include error checking. All given
%       %         transformation correspondences must be the same length
%       %         (i.e. numel(H_e2o) = numel(H_t2w) = numel(H_b2w)).
%
%       % Pre-process data
%       for i = 1:numel(H_t2w)
%           H_t2b{i} = invSE(H_b2w{i})*H_t2w{i};
%       end
%
%       % Solve for fixed transforms
%       [H_o2b,H_e2t,stats] = solveFixedTransforms(H_t2b,H_e2o);
%
%       % Calculate inverse transforms
%       H_b2o = invSE(H_o2b);
%       H_t2e = invSE(H_e2t);
%       % -----------------------------------------------------------------
%
%   M. Kutzer, 07Sep2022, USNA

% Update(s)
%   09Sep2022 - Updated to use parseVarargin_ZERO_fast
%   25Oct2022 - Updated documentation
%   27Oct2022 - Added surpressed output detection

%% Default options
ZERO = 1e-8;
fast = false;

%% Check input(s)
narginchk(2,4);

if ~iscell(H_x2a)
    error('H_x2a must be a cell array containing valid elements of SE(N).');
end

if ~iscell(H_y2b)
    error('H_y2b must be a cell array containing valid elements of SE(N).');
end

n_x2a = numel(H_x2a);
n_y2b = numel(H_y2b);
if n_x2a ~= n_y2b
    error('H_x2a and H_y2b must contain the same number of corresponding elements.')
end

% Define total number of correspondence pairs
n = n_x2a;

% Parse ZERO and "fast" values
[ZERO,fast,cellOut] = parseVarargin_ZERO_fast(varargin,ZERO,fast);

% TODO - check cellOut values for unused terms

%% Check for surpressed outputs using detectOutputSuppression.m
if exist('detectOutputSuppression.m','file') == 2
    % Function exists in the MATLAB path
    tfTilde = detectOutputSuppression(nargout);
else
    % Function does not exist in the MATLAB path
    tfTilde = false(1,nargout);
end
% Force 3-elment tfTilde
if numel(tfTilde) < 3
    tfTilde(3) = false;
end

%% Correct A/B pairs
str = {'H_x2a','H_y2b'};
if ~fast
    [H_x2a,H_y2b,info] = validCorrespondenceSE(H_x2a,H_y2b,ZERO);

    % Display removed pair information and altered transforms
    for j = 1:size(info.RemoveBin,2)
        % Display removed pair information
        if info.RemoveBin(1,j)
            fprintf('REMOVED PAIR: H_x2a{%d} - %s, H_y2b{%d} - %s\n',...
                j,info.RemoveMsg{1,j},j,info.RemoveMsg{2,j});
            continue
        end

        % Display altered transform information
        for i = 1:size(info.AlteredBin,1)
            if info.AlteredBin(i,j)
                fprintf('ALTERED TRANSFORM: %s{%d} - %s\n',...
                    str{i},j,info.RemoveMsg{i,j});
            end
        end
    end

    % Update number of pairs
    n = numel(H_x2a);
end

%% Redefine "fast" (we are assuming values are valid)
fastLocal = true;

%% Solve for H_b2a using AX = XB
% For an explanation of the method, see "EW450 Lectures\13 - Robot-Camera
% Calibration\Robot-Camera Calibration - Definition & Use.pptx"
%
%      A       X    =    X       B
%   H_ai2aj H_bi2ai = H_bj2aj H_bi2bj
%
%   Where H_bi2ai = H_bj2aj for all i,j

if nargout > 0 && ~tfTilde(1)
    % Initialize parameters
    iter = 0;
    nAB = (n^2 - n)/2;
    A = cell(1,nAB);
    B = cell(1,nAB);

    % Define i/j pairs
    for i = 1:n
        for j = 1:n
            % Isolate unique, non-identity transformation pairs
            %   We are keeping the upper-triangular portion of H_ai2aj and
            %   H_bi2bj to avoid values equal to the identity, and values that
            %   are the inverse of others.
            if i ~= j && i < j
                % H_ai2aj
                H_x2ai = H_x2a{i};
                H_x2aj = H_x2a{j};
                H_ai2aj = H_x2aj * invSE(H_x2ai,fastLocal);

                % H_bi2bj
                H_y2bi = H_y2b{i};
                H_y2bj = H_y2b{j};
                H_bi2bj = H_y2bj * invSE(H_y2bi,fastLocal);

                iter = iter+1;
                A{iter} = H_ai2aj;
                B{iter} = H_bi2bj;
            end
        end
    end
    fprintf('Number of A/B pairs: %d\n',iter);

    % Solve A * X = X * B
    X = solveAXeqXBinSE(A,B,ZERO,fastLocal);
    H_b2a = X;
    [tf,msg] = isSE(H_b2a,ZERO);
    if ~tf
        fprintf(...
            ['Value calculated for H_b2a is not a valid element of SE(3):\n\n',...
            '%s\nReplacing H_b2a with nearest element of SE(3).\n'],msg);
        H_b2a = nearestSE(H_b2a);
    end

    % Rename A/B matrices for approximating error
    H_ai2aj = A;
    H_bi2bj = B;
else
    % Output surpressed
    H_b2a = [];
end

%% Approximate H_b2a error
%   H_ai2aj H_bi2ai = H_bj2aj H_bi2bj
%       H_bi2aj     =     H_bi2aj

if nargout > 2 && ~tfTilde(1) && ~tfTilde(3)
    H_bi2ai = H_b2a;
    H_bj2aj = H_b2a;
    for i = 1:numel(H_ai2aj)
        LHS_H_bi2aj = H_ai2aj{i}*H_bi2ai;
        RHS_H_bi2aj = H_bj2aj*H_bi2bj{i};

        H_bi2bi{i} = invSE(LHS_H_bi2aj,fastLocal) * RHS_H_bi2aj;
        H_aj2aj{i} = LHS_H_bi2aj * invSE(RHS_H_bi2aj,fastLocal);
    end

    % Calculate mean
    % NOTE: For low error, these matrices should be very close to the identity
    muH_bi2bi = meanSE(H_bi2bi,ZERO,fastLocal);
    muH_aj2aj = meanSE(H_aj2aj,ZERO,fastLocal);
    % Calculate covariance
    % NOTE: For low error, these matrices should contain all values near zero
    SigmaH_bi2bi = covSE(H_bi2bi,muH_bi2bi,ZERO,fastLocal);
    SigmaH_aj2aj = covSE(H_aj2aj,muH_aj2aj,ZERO,fastLocal);
end

%% Solve for H_y2x using AX = XB
% For an explanation of the method, see "EW450 Lectures\13 - Robot-Camera
% Calibration\Robot-Camera Calibration - Definition & Use.pptx"
%
%      A       X    =    X       B
%   H_xi2xj H_yi2xi = H_yj2xj H_yi2yj
%
%   Where H_yi2xi = H_yj2xj for all i,j

if nargout > 1 && ~tfTilde(2)
    % Reinitialize parameters
    iter = 0;
    nAB = (n^2 - n)/2;
    A = cell(1,nAB);
    B = cell(1,nAB);

    % Define i/j pairs
    for i = 1:n
        for j = 1:n
            % Isolate unique, non-identity transformation pairs
            %   We are keeping the upper-triangular portion of H_ai2aj and
            %   H_bi2bj to avoid values equal to the identity, and values that
            %   are the inverse of others.
            if i ~= j && i < j
                % H_ai2aj
                H_xi2a = H_x2a{i};
                H_xj2a = H_x2a{j};
                H_xi2xj = invSE(H_xj2a,fastLocal) * H_xi2a;

                % H_bi2bj
                H_yi2b = H_y2b{i};
                H_yj2b = H_y2b{j};
                H_yi2yj = invSE(H_yj2b,fastLocal) * H_yi2b;

                iter = iter+1;
                A{iter} = H_xi2xj;
                B{iter} = H_yi2yj;
            end
        end
    end
    fprintf('Number of A/B pairs: %d\n',iter);

    % Solve A * X = X * B
    X = solveAXeqXBinSE(A,B,ZERO,fastLocal);
    H_y2x = X;
    [tf,msg] = isSE(H_y2x,ZERO);
    if ~tf
        fprintf(...
            ['Value calculated for H_b2a is not a valid element of SE(3):\n\n',...
            '%s\nReplacing H_b2a with nearest element of SE(3).\n'],msg);
        H_y2x = nearestSE(H_y2x);
    end

    % Rename A/B matrices for approximating error
    H_xi2xj = A;
    H_yi2yj = B;
else
    % Output surpressed
    H_y2x = [];
end

%% Approximate H_b2a error
%   H_xi2xj H_yi2xi = H_yj2xj H_yi2yj
%       H_yi2xj     =     H_yi2xj

if nargout > 2 && ~tfTilde(2) && ~tfTilde(3)
    H_yi2xi = H_y2x;
    H_yj2xj = H_y2x;
    for i = 1:numel(H_xi2xj)
        LHS_H_yi2xj = H_xi2xj{i}*H_yi2xi;
        RHS_H_yi2xj = H_yj2xj*H_yi2yj{i};

        H_yi2yi{i} = invSE(LHS_H_yi2xj,fastLocal) * RHS_H_yi2xj;
        H_xj2xj{i} = LHS_H_yi2xj * invSE(RHS_H_yi2xj,fastLocal);
    end

    % Calculate mean
    % NOTE: For low error, these matrices should be very close to the identity
    muH_yi2yi = meanSE(H_yi2yi,ZERO,fastLocal);
    muH_xj2xj = meanSE(H_xj2xj,ZERO,fastLocal);
    % Calculate covariance
    % NOTE: For low error, these matrices should contain all values near zero
    SigmaH_yi2yi = covSE(H_yi2yi,muH_yi2yi,ZERO,fastLocal);
    SigmaH_xj2xj = covSE(H_xj2xj,muH_xj2xj,ZERO,fastLocal);
end

%% Package output
if nargout > 2 && ~tfTilde(3)
    if ~tfTilde(1)
        stats.muH_b2b = muH_bi2bi;
        stats.muH_a2a = muH_aj2aj;
        stats.SigmaH_b2b = SigmaH_bi2bi;
        stats.SigmaH_a2a = SigmaH_aj2aj;
        stats.H_bi2bj = H_bi2bj;
        stats.H_ai2aj = H_ai2aj;
    end

    if ~tfTilde(2)
        stats.muH_y2y = muH_yi2yi;
        stats.muH_x2x = muH_xj2xj;
        stats.SigmaH_y2y = SigmaH_yi2yi;
        stats.SigmaH_x2x = SigmaH_xj2xj;
        stats.H_yi2yj = H_yi2yj;
        stats.H_xi2xj = H_xi2xj;
    end

    if tfTilde(1) && tfTilde(2)
        stats = [];
    end

    varargout{1} = stats;
end
