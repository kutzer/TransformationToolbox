function [tform,arg,tformStr,argStr] = parseKinematicString(kStr)
% PARSEKINEMATICSTRING parses a string or character array of rigid body
% transformations primitives (e.g., Rx, Ry, Rz, Tx, Ty, and Tz) and their
% associated magnitudes.
%
%   [tform,arg,tformStr,argStr] = parseKinematicString(kStr)
%
%   Input(s)
%       kStr - character array containing a "kinematic string" (ordered set
%              of rigid body motion primitizes and associated magnitudes). 
%
%           Example:
%               kStr = 'Rx(pi/2)Ty(5)Rz(0.3)' or
%               kStr = 'Rx(pi/2)*Ty(5)*Rz(0.3)
%
%   Output(s)
%       tform - 1xN cell array containing transformation function handles
%         arg - 1xN array containing transformation arguments
%    tformStr - 1xN cell array containing transformation functions defined 
%               as character arrays
%      argStr - 1xN cell array containing transformation arguments defined
%               as character arrays
%
%   Example:
%       kStr = 'Rz(pi/12)Tx((4-3)/18)Ty(52)Rx(-pi/4)Ty(23)';
%       [tform,arg,~,~] = parseKinematicString(kStr);
%       H_a2b = eye(4);
%       for i = 1:numel(tform)
%           H_a2b = H_a2b*tform{i}(arg(i));
%       end
%
%   Known Issues:
%       (1) This parsing code does not check for Shx vs Shxy arguments and
%           may result in errors when evaluating the resulting
%           transformations because of a dimension mismatch.
%
%   M. Kutzer, 14Jan2026, USNA

tfDebug = true;

%% Check input(s)
narginchk(1,1);
% TODO - check input(s)

%% Convert input to character array
kStr = char( kStr );

%% Remove '*' instances
kStr = strrep(kStr,'*','');

%% Parse transforms
% Regular Expression Breakdown:
% ([RTS][a-z]*)  -> Capture 1 or more letters starting with an R, T, or S (the transformation)
% \(           -> Match the literal opening parenthesis
% ( [^)]+ )    -> Capture everything that is NOT a closing parenthesis (the argument)
% \)           -> Match the literal closing parenthesis
%expression = '(?<op>[a-zA-Z]+)\((?<mag>[^)]+)\)';
expression = '(?<tforms>[RTS][a-z]*)\((?<arg>[^)]+)\)';

% Execute the match
matchStruct = regexp(kStr, expression, 'names');

% Isolate transforms
tformStr = {matchStruct.tforms};
n = numel(tformStr);

%% Parse transformation arguments
kWrk = kStr;

% Initialize output(s)
tform = cell(1,n);
arg = nan(1,n);
argStr = cell(1,n);

tfError = false;
for i = 1:n
    % Find argument indices
    i0 = strfind(kWrk,tformStr{i}  ) + numel(tformStr{i});
    if i < n
        i1 = strfind(kWrk,tformStr{i+1}) - 1;
    else
        i1 = numel(kWrk);
    end

    % Check transform instances
    if isempty(i0)
        % Bad string argument
        warning('%s( ) does not appear in the string.',tformStr{i});
        tfError = true;
        break
    end
    if isempty(i1)
        % Bad string argument
        warning('%s( ) does not appear in the string.',tformStr{i+1});
        tfError = true;
        break
    end

    % Isolate "first" instances
    if i < n
        if ~matches(tformStr{i},tformStr{i+1})
            i0 = i0(1);
            i1 = i1(1);
        else
            if numel(i1) < 2
                % Bad string argument
                warning('Second instance of %s( ) does not appear in the string.',tformStr{i+1});
                tfError = true;
                break
            end
            i0 = i0(1);
            i1 = i1(2);
        end
    end

    % Check for good argment
    if (i0-numel(tformStr{i})) ~= 1
        warning('Unexpected extra characters exist before transform: %s',kWrk( 1:(i0-numel(tformStr{i})) ));
        tfError = true;
        break
    end
    
    % Isolate argument
    kArg = kWrk(i0:i1);

    % Check argument
    i2 = strfind(kArg,')');
    if isempty(i2)
        % Bad string argument
        warning('%s( ) does not have a closing parenthesis.',tformStr{i});
        tfError = true;
        break
    end
    i2 = i2(end);
    
    % Check for good argment
    if i2 ~= numel(kArg)
        warning('Unexpected extra characters exist after argument: %s',kArg( (i2+1):end ));
        tfError = true;
        break
    end
    kArg = kArg(1:i2);
    
    % Append argument string
    argStr{i} = kArg;

    % Define numeric argument
    try
        arg(i) = eval(argStr{i});
    catch
        warning('%s%s is not valid',tformStr{i},kArg);
        tfError = true;
        break
    end
    
    % Define transformation function
    tform{i} = eval( sprintf('@(in)%s(in);',tformStr{i}) );

    if tfDebug
        fprintf('%s\n',kWrk);
        fprintf('%s - %s\n',tformStr{i},kArg);
    end

    % Remove parsed transformation
    kWrk = kWrk( (i1+1):end );
end

%% Check for error
if tfError
    tformStr = {};
    tform = {};
    arg = [];
    argStr = {};
end