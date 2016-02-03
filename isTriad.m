function [bin,kids] = isTriad(h)
%isTriad checks if the specified objects are valid "Triads".
%
%   bin - logical argument of "1" if specified object with the
%   corresponding index is a valid triad object and "0" otherwise.
%   kids - cell array containing the xyz axes line objects associated with 
%   the specified triad.
%
%   See also hgtransform triad showTriad
%
%   M. Kutzer 13May2015, USNA

bin = true(size(h));
kids = cell(size(h));
for i = 1:numel(h)
    if ~ishandle(h(i)) 
        bin(i) = false;
        continue
    end
    if ~strcmpi(get(h(i),'Type'),'hgtransform')
        bin(i) = false;
        continue
    end
    ckids = get(h(i),'Children');
    kids{i} = [];
    axs_tags = {'X-Axis','Y-Axis','Z-Axis'};
    for j = 1:numel(axs_tags)
        idx = find(~cellfun(@isempty, strfind(get(ckids,'Tag'),axs_tags{j})));
        for k = 1:numel(idx)
            if strcmpi( get(ckids(idx(k)),'Type'), 'Line' )
                kids{i}(end+1) = ckids(idx(k));
            end
        end
    end
    if numel(kids{i}) > 3
        warning('More than 3 xyz axes found for object.');
    end
end