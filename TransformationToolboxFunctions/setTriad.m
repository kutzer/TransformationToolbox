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


%% Get triad axes
[bin,all_kids] = isTriad(h);

%% Update properties
for i = 1:numel(bin)
    if ~bin(i)
        warning(sprintf('Index $d is not a valid triad.',i));
        continue
    end
    
    kids = all_kids{i};
    for j = 1:2:numel(varargin)
        switch lower(varargin{j})
            case 'linestyle'
                set(kids,varargin{j},varargin{j+1});
            case 'linewidth'
                set(kids,varargin{j},varargin{j+1});
            case 'matrix'
                set(h,varargin{j},varargin{j+1});
            case 'parent'
                set(h,varargin{j},varargin{j+1});
            case 'scale'
                s = varargin{j+1};
                if numel(s) == 1
                    s = repmat(s,1,3);
                end
                if numel(s) ~= 3
                    error('The scaling factor must be a singular value or a 3-element array.');
                end
                for j = 1:numel(kids)
                    xdata = get(kids(j),'XData');
                    ydata = get(kids(j),'YData');
                    zdata = get(kids(j),'ZData');
                    set(kids(j),'XData',xdata*s(1),'YData',ydata*s(2),'ZData',zdata*s(3));
                end
            case 'tag'
                set(h,varargin{j},varargin{j+1});
            case 'visible'
                set(kids,'Visible','off');
            otherwise
                % TODO - add check for properties in line or hgtransform, and
                % update property accordingly.
                warning(sprintf('Ignoring "%s," unexpected property.',varargin{j}));
        end
    end
end