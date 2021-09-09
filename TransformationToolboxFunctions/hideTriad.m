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
%   09Sep2021 - Updated to match showTriad

%% Hide labels
h = hideTriadLabels(h);

%% Hide axes
axs_tags = {'X-Axis','Y-Axis','Z-Axis'};
txt_tags = {'X-Label','Y-Label','Z-label'};
for i = 1:numel(h)
    % Show axes
    for j = 1:numel(axs_tags)
        kid = findobj('Parent',h(i),'Tag',axs_tags{j},'Type','line');
        set(kid,'Visible','off');
    end
    % Show labels
    for j = 1:numel(txt_tags)
        kid = findobj('Parent',h(i),'Tag',txt_tags{j},'Type','text');
        set(kid,'Visible','off');
    end
end