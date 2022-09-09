function H = smoothSE3(varargin)
% smoothSE3 Smooth cell array elements of SE3 sampled at a specific time
% interval.
%   H = smoothSE3(t,H) smooths SE3 data H using a 5-point moving average.
%
%   H = smoothSE3(t,H,SPAN) smooths SE3 data H using SPAN as the number of 
%   points used to compute each element of Z.
%  
%   H = smoothSE3(t,H,SPAN,METHOD) smooths SE3 data H with specified 
%   METHOD. The available methods are:
%  
%           'moving'   - Moving average (default)
%           'lowess'   - Lowess (linear fit)
%           'loess'    - Loess (quadratic fit)
%           'sgolay'   - Savitzky-Golay
%           'rlowess'  - Robust Lowess (linear fit)
%           'rloess'   - Robust Loess (quadratic fit)
%  
%   H = smoothSE3(t,H,METHOD) uses the default SPAN 5.
% 
%   See also smooth vee wedge
%
%   M. Kutzer 31Dec2014, USNA

% Updates
%   27Feb2016 - Replaced logm with logSO to compensate for negative
%               eigenvalue limitation of logm
%   09Sep2022 - Added validCorrespondenceSE to check for valid elements of
%               H

plotsOn = false;  % useful for debugging and visualizing output

%% Check for non-standard functions
%TODO - check for vee and wedge

%% Check inputs
if nargin < 2
    error('Not enough input arguments');
end
if nargin > 5
    error('Too many input arguments');
end
if nargin >= 2
    t = varargin{1};
    H = varargin{2};
    span = 10;          % larger span than typical "smooth.m"
    method = 'moving';  % same method as default in "smooth.m"
end
if nargin >= 3
    if ischar(varargin{3})
        method = varargin{3};
    else
        span = varargin{3};
    end
end
if nargin == 4
    %TODO - nothing checks for two methods and/or two spans
    if ischar(varargin{4})
        method = varargin{4};
    else
        span = varargin{4};
    end
end

%% Check for valid list of values
[H,info] = validCorrespondenceSE(H);

idxRmv = find(info.RemoveBin(1,:));
idxAlt = find(info.AlteredBin(1,:));
if isempty(idxRmv)
    msg = sprintf('The following values of H were removed because they are invalid:\n\t[ ');
    idxStr = sprintf('%d ',idxRmv);
    msg = sprintf('%s%s]\n\n\tSee validCorrespondenceSE for more info.',msg,idxStr);
    warning(msg);
end

if isempty(idxAlt)
    msg = sprintf('The following values of H were altered using nearestSE:\n\t[ ');
    idxStr = sprintf('%d ',idxAlt);
    msg = sprintf('%s%s]\n\n\tSee validCorrespondenceSE for more info.',msg,idxStr);
    warning(msg);
end

%% Check lengths
n = numel(H);
if numel(t) ~= n
    error('Each element of SE3 must correstpond to a specific time');
end

%% Retrieve data (SE3 -> R3 & SE3 -> so3)
v = zeros(3,n);
x = zeros(3,n);
for i = 1:n
    R = H{i}(1:3,1:3);
    v(:,i) = real( vee( logSO(R) ) );
    x(:,i) = real( H{i}(1:3,4) );
end

%% Smooth data
v_smooth = zeros(3,n);
x_smooth = zeros(3,n);
for i = 1:3
    v_smooth(i,:) = real( smooth(t,v(i,:),span,method)' );
    x_smooth(i,:) = real( smooth(t,x(i,:),span,method)' );
end

%% Plot results
if plotsOn
    fig = figure;
    for i = 1:6
        axs(i) = subplot(2,3,i,'Parent',fig);
        xlim([t(1),t(n)]);
        hold(axs(i),'on');
        if i < 4
            plot(axs(i),t,v(i,:),'.b');
            plot(axs(i),t,v_smooth(i,:),'m','LineWidth',1.5);
        else
            plot(axs(i),t,x(i-3,:),'.b');
            plot(axs(i),t,x_smooth(i-3,:),'m','LineWidth',1.5);
        end
    end
end

%% Package smooth data (R3 & so3 -> SE3)
for i = 1:n
    R = expm( wedge(v_smooth(:,i)) );
    H{i}(1:3,1:3) = R;
    H{i}(1:3,4) = x_smooth(:,i);
    H{i}(4,4) = 1;
end