function [h_j2p,fcn] = plotJoint(varargin)
% PLOTJOINT plots a specified joint type using common kinematic
% representations for prismatic and revolute joints.
%   [h_j2p,fcn] = plotJoint(jType,jAxis,jScale)
%
%   ___ = plotJoint(axs,___)
%
%   Input(s)
%       axs    - [OPTIONAL] handle defining the parent of the plotted joint
%       jType  - Joint type {'Revolute','Prismatic'}
%       jAxis  - 3x1 unit vector defining positive axis of translation or
%                rotation for the joint
%       jScale - positive scalar defining scale of joint
%
%   Output(s)
%       h_j2p - hgtransform object defining the "joint frame" relative to
%               the parent frame.
%       fcn   - anonymous function to rotate or translate the joint
%
%   M. Kutzer, 31Jan2025, USNA

%% Check input(s)
narginchk(3,4)

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

% Define type
if nargin >= nargin_i
    jType = lower( varargin{nargin_i} );
    nargin_i = nargin_i+1;
end

% Define axis
if nargin >= nargin_i
    jAxis = varargin{nargin_i};
    nargin_i = nargin_i+1;
end

% Define scale
if nargin >= nargin_i
    jScale = varargin{nargin_i};
    nargin_i = nargin_i+1;
end

% Check axis
if numel(jAxis) ~= 3
    error('Axis must be a 3x1 unit-vector.');
end
jAxis = reshape(jAxis,[],1);
jAxis = jAxis./norm(jAxis);

% Check scale
if numel(jScale) ~= 1
    error('Scale must be a positive scalar value.');
end
if jScale <= 0
    error('Scale must be a positive scalar value.');
end

% TODO - check remaining inputs

%% Create joint frame
h_j2p = triad('Parent',axs,'Scale',1.1*jScale,'LineWidth',1.5);

%% Plot joint
switch jType
    case 'revolute'
        % Define cylinder
        cyfit.Center = zeros(3,1);
        cyfit.Normal = jAxis;
        cyfit.Height = jScale;
        cyfit.Radius = 0.8*jScale/2;

        % Plot cylinder
        [ptc_cy,plt_cy] = plotCylinder(axs,cyfit);

        % Adjust plot object settings
        set(ptc_cy,'FaceColor','k','FaceAlpha',0.5,'EdgeColor','none');
        set(plt_cy,'Color','k','LineWidth',2);

        % Create function
        fcn = @(theta) revoluteJnt(h_j2p,jAxis,theta);
    case 'prismatic'
        % Define orientation
        R = vectorToSO(jAxis,3);

        % Define cube
        cufit.Center = zeros(3,1);
        cufit.Rotation = R;
        cufit.Dimension = jScale;

        % Plot cube
        ptc_cu = plotCube(axs,cufit);

        % Define square
        sqfit.Center = 1.05*jAxis*jScale;
        sqfit.Rotation = R;
        sqfit.Dimension = jScale;

        % Plot square
        ptc_sq = plotSquare(axs,sqfit);

        % Adjust plot object settings
        set(ptc_cu,'FaceColor','k','FaceAlpha',0.5,'EdgeColor','k','LineWidth',2);
        set(ptc_sq,'FaceColor','k','FaceAlpha',0.5,'EdgeColor','k','LineWidth',2);

        % Plot length
        plt_ln = plot(axs,nan,nan,'k','LineWidth',2);

        % Create function
        fcn = @(d) prismaticJnt(h_j2p,plt_ln,jAxis,jScale,d);
    otherwise
        error('Joint type must be revolute or prismatic.');
end

end

% Internal function(s)

function revoluteJnt(h_j2p,jAxis,theta)

set(h_j2p,'Matrix',expm(wedgeSE([jAxis*theta; zeros(3,1)])));

end

function prismaticJnt(h_j2p,plt_ln,jAxis,jScale,d)

set(h_j2p,'Matrix',Tx(jAxis(1)*d)*Ty(jAxis(2)*d)*Tz(jAxis(3)*d));
set(plt_ln,...
    'XData',[1.05*jAxis(1)*jScale, jAxis(1)*d],...
    'YData',[1.05*jAxis(2)*jScale, jAxis(2)*d],...
    'ZData',[1.05*jAxis(3)*jScale, jAxis(3)*d]);

end
