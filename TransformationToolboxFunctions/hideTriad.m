function h = hideTriad(h)
% HIDETRIAD hides the x, y, and z-direction and labels of the specified
% triad object (a transform object with three orthogonal lines as 
% children).
%   HIDETRIAD(h) hides the visualization of the axes associated with the
%   transform object(s) specified in h. Multiple objects must be specified
%   in an array.
%
%   See also hgtransform triad showTriad showTriadLabels hideTriadLabels
%
%   M. Kutzer 19Dec2014, USNA

%Updates
%   13May2015 - Updated definition for multiple triads and extended 
%               documentation
%   13Aug2015 - Updated to include triad labels

%% Hide labels
h = hideTriadLabels(h);

%% Hide axes
axs_tags = {'X-Axis','Y-Axis','Z-Axis'};
for i = 1:numel(h)
    kids = get(h(i),'Children');
    for j = 1:numel(axs_tags)
        idx = find(~cellfun(@isempty, strfind(get(kids,'Tag'),axs_tags{j})));
        for k = 1:numel(idx)
            if strcmpi( get(kids(idx(k)),'Type'), 'Line' )
                set(kids(idx(k)),'Visible','off');
            end
        end
    end
end