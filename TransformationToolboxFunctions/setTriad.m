function setTriad(h,varargin)
%setTriad sets the specified properties of one or more triad objects
%
%   Properties:
%       'linestyle' - [ {-} | -- | : | -. | none ]
%       'linewidth' - Used to define thickness lines representing axes
%                     Default value is 0.50
%                     Value must be finite and greater than zero
%       'matrix'    - 4x4 homogenious transform
%       'parent'    - Axes or other hgtransform handle
%       'scale'     - Used to define the length of each axes line
%                     Default value is 1.00
%                     Value(s) must be finite and greater than zero
%                     A scalar value scales the x, y, and z-axis equally
%                     A 3-element array (e.g. [1,2,3]) scales each of the
%                       axes seperately.
%       'tag'       - String describing object
%       'visible'   - [ {on} | off ]
%                     This property hides the visualization, not the
%                       hgtransform object.
%
%   See also triad isTriad getTriad showTriad hideTriad
%
%   M. Kutzer 19Dec2014, USNA

% Updates:
%   13Jul2020 - Corrected scaling update.
%   25Apr2022 - Added AxisLabels

%% Get triad axes
[bin,kAxesALL,kLabelsALL] = isTriad(h);

%% Update properties
for i = 1:numel(bin)
    if ~bin(i)
        warning(sprintf('Index $d is not a valid triad.',i));
        continue
    end

    kAxes = kAxesALL{i};
    kLabels = kLabelsALL{i};
    for j = 1:2:numel(varargin)
        switch lower(varargin{j})
            case 'linestyle'
                set(kAxes,varargin{j},varargin{j+1});
            case 'linewidth'
                set(kAxes,varargin{j},varargin{j+1});
            case 'matrix'
                set(h(i),varargin{j},varargin{j+1});
            case 'parent'
                set(h(i),varargin{j},varargin{j+1});
            case 'scale'
                s = varargin{j+1};
                if numel(s) == 1
                    s = repmat(s,1,3);
                end
                if numel(s) ~= 3
                    error('The scaling factor must be a singular value or a 3-element array.');
                end
                % Update axes
                for k = 1:numel(kAxes)
                    X(1,:) = get(kAxes(k),'XData');
                    X(2,:) = get(kAxes(k),'YData');
                    X(3,:) = get(kAxes(k),'ZData');
                    dX = abs( diff(X,1,2) );
                    idx = find(dX == max(dX),1,'first');
                    switch idx
                        case 1
                            set(kAxes(k),...
                                'XData',[0,1]*s(1),...
                                'YData',[0,0],...
                                'ZData',[0,0]);
                        case 2
                            set(kAxes(k),...
                                'XData',[0,0],...
                                'YData',[0,1]*s(2),...
                                'ZData',[0,0]);
                        case 3
                            set(kAxes(k),...
                                'XData',[0,0],...
                                'YData',[0,0],...
                                'ZData',[0,1]*s(3));
                        otherwise
                            error('Something went wrong!');
                    end
                end
                % Update labels
                for k = 1:numel(kLabels)
                    kPos = abs( get(kLabels(k),'Position') );
                    idx = find(kPos == max(kPos),1,'first');
                    switch idx
                        case 1
                            set(kLabels(k),'Position',[1,0,0]*s(1));
                        case 2
                            set(kLabels(k),'Position',[0,1,0]*s(2));
                        case 3
                            set(kLabels(k),'Position',[0,0,1]*s(3));
                        otherwise
                            error('Something went wrong!');
                    end
                end
            case 'tag'
                set(h,varargin{j},varargin{j+1});
            case 'visible'
                set(kAxes,'Visible','off');
            case 'axislabels'
                axislabels = varargin{j+1};
                if numel(axislabels) ~= numel(kAxes)
                    error('An axis label must be specified for each axis.')
                end
                if ~iscell(axislabels)
                    error('Axis labels must be specified as a 1x3 cell array.');
                end

                if isempty(kLabels)
                    for k = 1:numel(kAxes)
                        txt(k) = text('Parent',h(i));
                    end
                else
                    txt = kLabels;
                end

                axislabel_tags = {'X-Label','Y-Label','Z-Label'};
                for k = 1:numel(axislabels)
                    % Convert label to string argument
                    if ~ischar(axislabels{k})
                        axislabels{k} = num2str(axislabels{k});
                    end
                    % Get label position
                    xdata = get(kAxes(k),'XData');
                    ydata = get(kAxes(k),'YData');
                    zdata = get(kAxes(k),'ZData');
                    % Set label position
                    set(txt(k),'Position',...
                        [xdata(end),ydata(end),zdata(end)]);
                    % Set label properties
                    %{
                    set(txt(k),'Parent',h(i),...
                        'HorizontalAlignment','Left',...
                        'VerticalAlignment','Bottom',...
                        'Color',[0 0 0],...
                        'FontName','Helvetica',...
                        'FontSize',10,...
                        'Tag',axislabel_tags{k});
                    %}
                    set(txt(k),'Parent',h(i),...
                        'HorizontalAlignment','Left',...
                        'VerticalAlignment','Bottom',...
                        'Tag',axislabel_tags{k},...
                        'String',axislabels{k});
                end
            otherwise
                % kAxesALL,kLabelsALL
                if isprop(h(i),varargin{j})
                    set(h(i),varargin{j},varargin{j+1});
                elseif all(isprop(kAxes,varargin{j}),'all')
                    set(kAxes,varargin{j},varargin{j+1});
                elseif all(isprop(kLabels,varargin{j}),'all')
                    set(kLabels,varargin{j},varargin{j+1});
                else
                    % TODO - add check for properties in line or hgtransform, and
                    % update property accordingly.
                    warning(sprintf('Ignoring "%s," unexpected property.',varargin{j}));
                end
        end
    end
end