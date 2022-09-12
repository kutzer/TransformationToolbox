function [ZERO,fast,cellOut] = parseVarargin_ZERO_fast(cellIn,ZERO,fast)
% PARSEcellIn_ZERO_FAST parses a variable input argument for ZERO and
% fast parameters.
%   [ZERO,fast,cellOut] = parsecellIn_ZERO_fast(cellIn,ZERO,fast)
%
%   Input(s)
%       cellIn - n-element cell array that may include values for ZERO
%                and/or fast parameters. Additional parameters can be
%                included and will be returned in cellOut.
%       ZERO   - [OPTIONAL] positive scalar *default* value that is
%                sufficiently close to zero to be assumed zero
%                (e.g. ZERO = 1e-8). If a default "ZERO" is not specified,
%                a default of ZERO = [] is used.
%       fast   - [OPTIONAL] true/false logical *default* value indicating
%                whether to skip checking a specified property or
%                properties. If "fast" is not specified, a default of
%                fast = false is used.
%           fast = true    - Skip checking property or properties
%           fast = [false] - Check property or properties
%
%   Output(s)
%       ZERO    - positive scalar value that is sufficiently close to zero
%                 to be assumed zero (e.g. ZERO = 1e-8). If "ZERO" is not
%                 specified, the default value is used.
%       fast    - true/false logical value indicating whether to
%                 skip checking a specified property or properties. If
%                 "fast" is not specified, the default value is used.
%           fast = true  - Skip checking property or properties
%           fast = false - Check property or properties
%       cellOut - remaining values contained in cellIn that do not meet the
%                 ZERO or fast criteria
%
%   NOTE: This function is used primarily in the Transformation Toolbox for
%         parsing optional, variable input arguments that are used across a
%         wide range of functions.
%
%   M. Kutzer, 09Sep2022, USNA

% Update(s)
%   12Sep2022 - Accounted for no inputs

%% Set default values
if nargin < 3
    fast = false;
end
if nargin < 2
    ZERO = [];
end
cellOut = {};

%% Check number of inputs
narginchk(1,3);

%% Check cellIn
if ~iscell(cellIn)
    error('Input must be a cell array');
end

if isempty(cellIn)
    % No values provided
    return
end

if iscell(cellIn{1})
    % varargin as input workaround
    cellIn = cellIn{1};
end

%% Parse cellIn
isDefined_ZERO = false; % Flag noting whether a value for "ZERO" is found
isDefined_fast = false; % Flag noting whether a value for "fast" is found
tfRemove = false(size(cellIn));
for i = 1:numel(cellIn)
    %fprintf('Input %d: ',i);
    val_i = cellIn{i};

    % NUMERIC VALUE(S)
    if isnumeric(val_i)
        switch numel(val_i)
            case 0
                %fprintf('1\n');
                % Candidate "ZERO" is specified as []

                if isDefined_ZERO
                    if isempty(ZERO)
                        warning('Multiple candidate values for the "ZERO" parameter are present. Using first value: "ZERO = []".');
                    else
                        warning('Multiple candidate values for the "ZERO" parameter are present. Using first value: "ZERO = %f".',ZERO);
                    end
                else
                    ZERO = val_i;
                    isDefined_ZERO = true;
                    tfRemove(i)    = true;
                end
                continue
            case 1
                % Possible "ZERO" or "fast" value

                if val_i == 0 || val_i == 1
                    %fprintf('2\n');
                    % Candidate "fast" value specified as 0 or 1 [LEGACY]

                    % TODO - Warn user that "fast" will be logical values
                    %        only in the future
                    if isDefined_fast
                        warning('Multiple candidate values for the "fast" parameter are present. Using first value: "fast = logical(%d)".',fast);
                    else
                        fast = logical(val_i);
                        isDefined_fast = true;
                        tfRemove(i)    = true;
                    end
                    continue
                end

                if val_i >= 0
                    %fprintf('3\n');
                    % Candidate "ZERO" value is specified
                    if isDefined_ZERO
                        if isempty(ZERO)
                            warning('Multiple candidate values for the "ZERO" parameter are present. Using first value: "ZERO = []".');
                        else
                            warning('Multiple candidate values for the "ZERO" parameter are present. Using first value: "ZERO = %f".',ZERO);
                        end
                    else
                        ZERO = val_i;
                        isDefined_ZERO = true;
                        tfRemove(i)    = true;
                    end
                    continue
                end
        end
        continue
    end

    % LOGICAL, CHARACTER, & STRING VALUE(S)
    switch lower( class(val_i) )
        case 'logical'
            if numel(val_i) == 1
                %fprintf('4\n');
                % Candidate "fast" value is specified as logical scaler
                if isDefined_fast
                    warning('Multiple candidate values for the "fast" parameter are present. Using first value: "fast = logical(%d)".',fast);
                else
                    fast = val_i;
                    isDefined_fast = true;
                    tfRemove(i)    = true;
                end
                continue
            end
        case {'char','string'}
            switch lower( char(val_i) )
                case 'fast'
                    %fprintf('5\n');
                    % Candidate "fast" value specified as character
                    % array [LEGACY]

                    % TODO - Warn user that "fast" will be logical
                    %        values only in the future
                    if isDefined_fast
                        warning('Multiple candidate values for the "fast" parameter are present. Using first value: "fast = logical(%d)".',fast);
                    else
                        fast = true;
                        isDefined_fast = true;
                        tfRemove(i)    = true;
                    end
                    continue
            end
    end
    %fprintf('\n')
end

%% Check ZERO and fast values
if ~isDefined_ZERO
    if ~isnumeric(ZERO) || numel(ZERO) > 1
        error('"ZERO" must be defined as a scalar value greater than 0 or [].');
    end
    if numel(ZERO) == 1
        if ZERO(1) < 0
            error('"ZERO" must be defined as a scalar value greater than 0 or [].');
        end
    end
end

if ~isDefined_fast
    % Account for legacy cases of user supplied "fast" values
    switch lower( class(fast) )
        case {'char','string'}
            switch lower( char(fast) )
                case 'fast'
                    % Candidate "fast" value specified as character
                    % array [LEGACY]

                    % TODO - Warn user that "fast" will be logical
                    %        values only in the future
                    fast = true;
            end
        otherwise
            if isnumeric(fast) || numel(fast) == 1
                % Candidate "fast" value specified as non-logical
                % scalar [LEGACY]

                % TODO - Warn user that "fast" will be logical
                %        values only in the future
                if fast == 0 || fast == 1
                    fast = logical(fast);
                end
            end
    end

    if ~islogical(fast) || numel(fast) ~= 1
        error('"fast" must be defined as logical scalar value.')
    end
end

%% Pass remaining non-parsed values
cellOut = cellIn(~tfRemove);