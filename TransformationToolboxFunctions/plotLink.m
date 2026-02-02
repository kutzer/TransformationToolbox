function varargout = plotLink(varargin)
% PLOTLINK plots a rigid link along a designated direction
%   plt = plotLink(link)
%
%   [h_l2p,plt] = plotLink(link)
%
%   plt = plotLink(link,offsets)
%
%   [h_l2p,plt] = plotLink(link,offsets)
%
%   ___ = plotLink(axs,___)
%
%   Input(s)
%       axs    - [OPTIONAL] handle defining the parent of the plotted link
%       link   - 3x1 array specifying link magnitude and direction
%       offset - [OPTIONAL] 1x2 array specifying the initial and final 
%                offsets for the link. Default value is offset = [0,0];
%
%   Output(s)
%       h_l2p - hgtransform object defining the "link frame" relative to
%               the parent frame.
%       plt   - line object visualizing link
%
%   M. Kutzer, 31Jan2025, USNA

% Updates
%   02Feb2026 - Added hgtransform output and corrected documentation

%% Check input(s)
narginchk(1,3)

% Define parent
parentGiven = false;
nargin_i = 1;
if ishandle(varargin{nargin_i})
    switch get(varargin{nargin_i},'Type')
        case 'axes'
            parentGiven = true;
        case 'hgtransform'
            parentGiven = true;
    end
end

if parentGiven
    axs = varargin{nargin_i};
    nargin_i = nargin_i+1;
else
    axs = gca;
end

switch get(axs,'Type')
    case 'axes'
        hold(axs,'on');
end

% Define link
if nargin >= nargin_i
    link = lower( varargin{nargin_i} );
    nargin_i = nargin_i+1;
end

% Define offset
offset = [0,0];
if nargin >= nargin_i
    offset = varargin{nargin_i};
    nargin_i = nargin_i+1;
end

% Check axis
if numel(link) ~= 3
    error('Link must be a 3x1 array.');
end
link = reshape(link,[],1);
v_mag = norm(link);
v_hat = link./v_mag;

% Check offset
if numel(offset) ~= 2
    error('Scale must be a 2-element vector.');
end

%% Plot 
X = [offset(1)*v_hat, (v_mag-offset(2))*v_hat];
plt = plot3(X(1,:),X(2,:),X(3,:),'-k','LineWidth',2,'Parent',axs);

if nargout <= 1
    varargout{1} = plt;
    return
end

%% Create frame
h_l2p = hgtransform('Parent',axs,'Matrix',Tx(link(1))*Ty(link(2))*Tz(link(3)));

varargout{1} = h_l2p;
varargout{2} = plt;

% TODO - check remaining inputs