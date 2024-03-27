function prop = getTriad(h,varargin)
% GETTRIAD gets properties specified as name/value pairs from an
% hgtransform object associated with a triad.
%
%   prop = getTriad(h,NAME)
%
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
%   M. Kutzer, 27Mar2024, USNA

%% Check input(s)
if nargin < 2
    error('Inputs must be an hgtransform object and name/value pairs');
end

%% Recover lines and text of triad
txt = findobj(h,'Parent',h,'Type','Text');
plt = findobj(h,'Parent',h,'Type','Line');
tag_plt = get(plt,'Tag');
tag_txt = get(txt,'Tag');

lbls = {'X-Axis','Y-Axis','Z-Axis'};

tf_plt = false(1,3);
idx_plt = zeros(1,3);
tf_txt = false(1,3);
idx_txt = zeros(1,3);
for i = 1:numel(lbls)
    tfMatch_plt = matches(tag_plt,lbls{i});
    switch nnz(tfMatch_plt)
        case 0
            % No axis line is associated
            error('No "%s" line object is associated with this hgtransform.',...
                lbls{i});
        case 1
            % One axis is associated
            tf_plt(i) = true;
            idx_plt(i) = find(tfMatch_plt);
        otherwise
            % Multiple axis lines are associated
            error('Multiple "%s" line objects are associated with this hgtransform.',...
                lbls{i});
    end

    tfMatch_txt = matches(tag_txt,lbls{i});
    switch nnz(tfMatch_txt)
        case 0
            % No axis line is associated
            error('No "%s" text object is associated with this hgtransform.',...
                lbls{i});
        case 1
            % One axis is associated
            tf_txt(i) = true;
            idx_txt(i) = find(tfMatch_txt);
        otherwise
            % Multiple axis lines are associated
            error('Multiple "%s" text objects are associated with this hgtransform.',...
                lbls{i});
    end

end

if all(tf_plt,'all')
    % Re-order line and text objects to x/y/z
    plt = plt(idx_plt);
    txt = txt(idx_txt);

    color = get(plt,'Color');
    str = get(txt,'String');
    for i = 1:numel(plt)
        X(1,:) = get(plt(i),'XData');
        X(2,:) = get(plt(i),'YData');
        X(3,:) = get(plt(i),'ZData');

        % TODO - check for 2 points
        sc(i) = norm( diff(X,1,2) );
    end
end

%% Recover properties
for i = 1:numel(varargin)
    if ~ischar(varargin{i}) && ~isstring(varargin{i})
        % Property is not a valid character array or string
        warning('Property %d is not a valid character array or string.');
        continue
    end

    if isstring(varargin{i})
        varargin{i} = convertStringsToChars(varargin{i});
    end

    switch lower(varargin{i})
        case 'line'
