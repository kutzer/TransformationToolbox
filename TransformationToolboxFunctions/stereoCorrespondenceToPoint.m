function [P_ca,P_cb,err_ma,err_mb] = stereoCorrespondenceToPoint(P_ma,P_mb,A_ca2ma,A_cb2mb,H_ca2cb,sb_estimate)
% STEREOCORRESPONDENCETOPOINT Estimate the 3D position of a point given a
% pixel correspondence from two cameras with known extrinsics and
% extrinsics. 
%   ___ = STEREOCORRESPONDENCETOPOINT(P_ma,P_mb,A_ca2ma,A_cb2mb,H_ca2cb)
%
%   ___ = STEREOCORRESPONDENCETOPOINT(P_ma,P_mb,A_ca2ma,A_cb2mb,H_ca2cb,sb)
%   estimates the 3D position of a point given:
%
%   Input(s)
%       P_ma    - 2D pixel coordinate of the point relative to camera a
%       P_mb    - 2D pixel coordinate of the point relative to camera b
%       A_ca2ma - Camera a intrisic matrix 
%       A_cb2mb - Camera b intrisic matrix 
%       H_ca2cb - Extrinsics defining camera a relative to camera b
%       sb      - [OPTIONAL] Initial depth estimate of the point relative 
%                 to camera b
%
%   [P_ca,P_cb] = STEREOCORRESPONDENCETOPOINT(___) returns
%
%   [P_ca,P_cb,err_ma,err_mb] = STEREOCORRESPONDENCETOPOINT(___) returns
%
%   Output(s)
%       P_ca    - The point defined relative to camera a
%       P_cb    - The point defined relative to camera b
%       err_ma  - The reprojection error of the P_ca estimate
%       err_mb  - The reprojection error of the P_cb estimate
%
%   M. Kutzer, 17Oct2019, USNA

% Updates
%   26Feb2021 - Updated documentation

%% Check inputs
% Number of inputs
narginchk(5,6);

% Pixel coordinate dimensions
if numel(P_ma) < 2 || numel(P_ma) > 3
    error('Pixel coordinate for camera a must be specified as a 2 or 3 element array.');
end

if numel(P_mb) < 2 || numel(P_mb) > 3
    error('Pixel coordinate for camera b must be specified as a 2 or 3 element array.');
end

% Intrisic dimensions
if ~ismatrix(A_ca2ma) || size(A_ca2ma,1) ~= 3 || size(A_ca2ma,2) ~= 3
    error('Intrinsic matrix for camera a must be 3x3.');
end

if ~ismatrix(A_cb2mb) || size(A_cb2mb,1) ~= 3 || size(A_cb2mb,2) ~= 3
    error('Intrinsic matrix for camera a must be 3x3.');
end

% Intrinsic matrix form
if A_ca2ma(3,1) ~= 0 || A_ca2ma(3,2) ~= 0 || A_ca2ma(3,3) ~= 1
    error(['Intrinsic matrix must be of the form\n',...
        '\t[ a11  a12  a13 ]\n',...
        '\t[ a21  a22  a23 ]\n',...
        '\t[   0    0    1 ]']);
end

if A_cb2mb(3,1) ~= 0 || A_cb2mb(3,2) ~= 0 || A_cb2mb(3,3) ~= 1
    error(['Intrinsic matrix must be of the form\n',...
        '\t[ a11  a12  a13 ]\n',...
        '\t[ a21  a22  a23 ]\n',...
        '\t[   0    0    1 ]']);
end

% Extrinsics
if ~isSE(H_ca2cb)
    error('Extrinsic matrix descibing the pose of camera a relative to camera b must be a valid 3D rigid body transform.');
end

%% Set default value of sb_estimate
if nargin < 6
    sb_estimate = 10*norm(H_ca2cb(1:3,4));
    warning('No estimated depth of the point relative to camera b was provided. Estimating depth as %.3f.',sb_estimate);
end

if numel(sb_estimate) ~= 1
    error('Estimated depth of the point relative to camera b must be specified as a scalar.');
end

%% Reshape the pixel coordinates and ensure they are homogeneous
% Make sure coordinates are column vectors
P_ma = reshape(P_ma,[],1);
P_mb = reshape(P_mb,[],1);
% Make coordinates homogeneous
P_ma(3,1) = 1;
P_mb(3,1) = 1;

%% Define debug flags
plotsOn = false;

if plotsOn
    fig_s = figure('Name','stereoCorrespondenceToPoint.m');
    axs_s = axes('Parent',fig_s);
    hold(axs_s,'on');
    title(axs_s,'Depth Convergence');
    xlabel(axs_s,'s_a');
    ylabel(axs_s,'s_b');
    
    ptc_sTol = patch('Parent',axs_s,'FaceColor','r','FaceAlpha',0.3,'EdgeColor','k');
    plt_sAll = plot(axs_s,0,0,'b','LineWidth',1.5);
    plt_sAvg = plot(axs_s,0,0,'*b');
end

%% Make sure pixel coordinates are column vectors
P_ma = reshape(P_ma,[],1);
P_mb = reshape(P_mb,[],1);

%% Convert pixel coordinates to scaled camera coordinates
% Define matrix frames relative to camera frames
%   $A_{m_a}^{c_a} = \left( A_{c_a}^{m_a} \right)^{-1}$
%   $A_{m_b}^{c_b} = \left( A_{c_b}^{m_b} \right)^{-1}$
A_ma2ca = A_ca2ma^(-1);
A_mb2cb = A_cb2mb^(-1);
% Define "scaled" points relative to camera frames
%   $\frac{1}{s_a} \vec{p}^{c_a} = A_{m_a}^{c_a} \vec{p}^{m_a}$
%   $\frac{1}{s_b} \vec{p}^{c_b} = A_{m_b}^{c_b} \vec{p}^{m_b}$
P_csa = A_ma2ca * P_ma;
P_csb = A_mb2cb * P_mb;

%% Recover scale (i.e. depth relative to each camera)
% TODO - Reference or prove this method

% Define camera b relative to camera a
%   $H_{c_b}^{c_a} = \left( H_{c_a}^{c_b} \right)^{-1}$
H_cb2ca = H_ca2cb^(-1);

% Guess scale/depth values relative to camera b
%   Assume 100 equally spaced values starting at 10% of the estimate and 
%   ending at 1000% of the estimate.
%   TODO - Validate this decision mathematically
sb = linspace(0.1,10.0,100).*sb_estimate;

% Define acceptable convergence value
convTolerance = 0.001;
% Initialze current convergence value
convValue = inf;
while convValue > convTolerance
    % Estimate 3D point relative to camera a using estimated depth to camera b
    %   $\vec{p}^{c_a} = s_b H_{c_b}^{c_a} \left( \begin{array}{c} \frac{1}{s_b} \vec{p}^{c_b} \\ 1 \end{array} \right)$
    %   $\vec{p}^{c_a} = s_b R_{c_b}^{c_a} \left( \frac{1}{s_b} \vec{p}^{c_b} \right) + \vec{d}_{c_b}^{c_a}$ 
    P_ca = sb.*repmat(H_cb2ca(1:3,1:3)*P_csb,1,numel(sb)) + repmat(H_cb2ca(1:3,4),1,numel(sb));
    % Extract estimates of depth to camera a
    sa = P_ca(3,:);
    
    % Estimate 3D point relative to camera b using estimated depth to camera a
    %   $\vec{p}^{c_b} = s_a H_{c_a}^{c_b} \left( \begin{array}{c} \frac{1}{s_a} \vec{p}^{c_a} \\ 1 \end{array} \right)$
    %   $\vec{p}^{c_b} = s_a R_{c_a}^{c_b} \left( \frac{1}{s_a} \vec{p}^{c_a} \right) + \vec{d}_{c_a}^{c_b}$ 
    P_cb = sa.*repmat(H_ca2cb(1:3,1:3)*P_csa,1,numel(sa)) + repmat(H_ca2cb(1:3,4),1,numel(sa));
    % Extract estimates of depth to camera b
    sb = P_cb(3,:);
    
    % Calculate convergence value
    %   TODO - Validate this decision mathematically
    convValue = max( [max(sa) - min(sa), max(sb) - min(sb)] );
    
    % Visualize convergence
    if plotsOn
        % Update depth to camera b vs depth to camera a
        set(plt_sAll,'XData',sa,'YData',sb);
        % Update average
        set(plt_sAvg,'XData',mean(sa),'YData',mean(sb));
        % Update tolerance box
        verts = convTolerance/2.*[-1,-1; 1, -1; 1, 1; -1, 1];
        verts = verts + ...
            [repmat( (max(sa) - min(sa))/2 + min(sa) ,4,1),...
             repmat( (max(sb) - min(sb))/2 + min(sb) ,4,1)];
        ptc_sTol.Vertices = verts;
        ptc_sTol.Faces = 1:4;
        % Update plot
        drawnow
    end
end
% Average the result
sa = mean(sa);
sb = mean(sb);

%% Calculate 3D position of point
% Relative to camera a
P_ca = sb.*repmat(H_cb2ca(1:3,1:3)*P_csb, 1,numel(sb)) + repmat(H_cb2ca(1:3,4),1,numel(sb));
% Relative to camera b
P_cb = sa.*repmat(H_ca2cb(1:3,1:3)*P_csa, 1,numel(sa)) + repmat(H_ca2cb(1:3,4),1,numel(sa));

% Return if the user is not requesting reprojection error
if nargout < 3 && ~plotsOn
    return
end

%% Reproject points back into camera space 
% Reproject point to camera a matrix frame 
P_ma_est = A_ca2ma * P_ca;
P_ma_est = P_ma_est./P_ma_est(3);
% Calculate camera a reprojection error
err_ma = norm(P_ma_est - P_ma);

% Reproject point to camera b matrix frame 
P_mb_est = A_cb2mb * P_cb;
P_mb_est = P_mb_est./P_mb_est(3);
% Calculate camera b reprojection error
err_mb = norm(P_mb_est - P_mb);

if plotsOn
    fig_m = figure('Name','stereoCorrespondenceToPoint.m');
    for i = 1:2
        axs_m(i) = subplot(1,2,i,'Parent',fig_m,'YDir','reverse');
        hold(axs_m(i),'on');
        daspect(axs_m(i),[1 1 1]);
        xlabel(axs_m(i),'x (pixels)');
        ylabel(axs_m(i),'y (pixels)');
    end
    title(axs_m(1),'Camera a');
    title(axs_m(2),'Camera b');
    
    plot(axs_m(1),P_ma(1),P_ma(2),'xr');
    plot(axs_m(2),P_mb(1),P_mb(2),'xb');
    plot(axs_m(1),P_ma_est(1),P_ma_est(2),'or');
    plot(axs_m(2),P_mb_est(1),P_mb_est(2),'ob');
end
