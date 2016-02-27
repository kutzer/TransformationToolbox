function DH_obj = plotDHtable(varargin)
% PLOTDHTABLE plots the coordinate frames and links associated with a DH
% table.
%
%       DHtable = [theta_0,d_0,a_0,alpha_0; 
%                  theta_1,d_1,a_1,alpha_1;
%                  ...
%                  theta_N,d_N,a_N,alpha_N];
%       H_i = Rz(theta_i)*Tz(d_i)*Tx(a_i)*Rx(alpha_i);
%       H = H_0*H_1*...*H_N;
%
%   See also DHtableToFkin DH
%   
%   M. Kutzer 16July2015, USNA

axs = [];
linkLbls = true;
if nargin > 0
    DHtable = varargin{1};
end
if nargin > 1
    axs = varargin{1};
    DHtable = varargin{2};
end
if nargin > 2
    props = varargin(3:end);
    for i = 1:2:numel(props)
        switch lower(props{i})
            case 'linklabels'
                switch lower(props{i+1})
                    case 'on'
                        linkLbls = true;
                    case 'off'
                        linkLbls = false;
                    otherwise
                        error('Unexpected property value.');
                end
            otherwise
                error('Unexpected property name.');
        end
    end
end

if isempty(axs)
    fig = figure('Name','DH Table Plot');
    axs = axes('Parent',fig);
    daspect([1 1 1]);
    hold on
end

%TODO - auto select the scale for triads given lengths of d and a
Scale = 50;
LineWidth = 2;

h(1) = triad('Parent',axs,'Scale',Scale,'LineWidth',LineWidth,...
    'Matrix',eye(4),'AxisLabels',{'x_0','y_0','z_0'});
N = size(DHtable,1);
for i = 1:N
    theta_i = DHtable(i,1);
    d_i     = DHtable(i,2);
    a_i     = DHtable(i,3);
    alpha_i = DHtable(i,4);
    H_i{i} = DH(theta_i,d_i,a_i,alpha_i);
	
    AxisLabels{1} = sprintf('x_%d',i);
    AxisLabels{2} = sprintf('y_%d',i);
    AxisLabels{3} = sprintf('z_%d',i);
    
    h(i+1) = triad('Parent',h(i),'Scale',Scale,'LineWidth',LineWidth,...
        'Matrix',H_i{i},'AxisLabels',AxisLabels);
    g = hgtransform('Parent',h(i+1),'Matrix',Rx(-alpha_i));
    plt(1) = plot3([0,-a_i],[0,0],[0,0],':k');
    plt(2) = plot3([-a_i,-a_i],[0,0],[0,-d_i],':k');
    set(plt,'Parent',g);
    
    if linkLbls
        if a_i ~= 0
            try
                txt = text(-(a_i/2),0,0,sprintf('a_%d',i));
                set(txt,'Parent',g);
            catch
                warning('Unable to display a_%d.',i);
            end
        end
        if d_i ~= 0
            try
                txt = text(-a_i,0,-d_i/2,sprintf('d_%d',i));
                set(txt,'Parent',g);
            catch
                warning('Unable to display d_%d.',i);
            end
        end
    end
    
end

view(3);
DH_obj = h;