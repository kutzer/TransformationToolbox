function vars = readFileToArray(fname)
% READFILETOARRAY reads a file created using "writeArrayToFile"
%   out = readFileToArray(fname)
%
%   Input(s)
%       fname - character array specifying filename containing array
%               information
%
%   Output(s)
%       vars - structured array containing fields matching variable names
%              specified in file
%
%   Example
%       fname = 'EXAMPLE_writeArrayToFile.txt';
%       vars = readFileToArray(fname);
%
%       'EXAMPLE_writeArrayToFile.txt' File Contents:
%           H_a2b{1} = reshape([___, ___, ___, ___],2,2);\n
%           H_b2c{1} = reshape([___, ___, ___, ___],2,2);\n
%           H_a2b{2} = reshape([___, ___, ___, ___],2,2);\n
%           H_b2c{2} = reshape([___, ___, ___, ___],2,2);\n
%           H_a2b{3} = reshape([___, ___, ___, ___],2,2);\n
%           H_b2c{3} = reshape([___, ___, ___, ___],2,2);\n
%
%       Description of "vars" output
%           vars.H_a2b - 3-element cell array containing 2x2 arrays
%           vars.H_b2c - 3-element cell array containing 2x2 arrays
%
%   See also writeArrayToFile
%
%   M. Kutzer, 14Oct2022, USNA

vars = [];

%% Check input(s)
switch class(fname)
    case {'char','string'}
        % Open file
        fID = fopen(fname,'r');
        if fID < 0 
            error('Unable to open "%s". Try a new filename.',fname);
        end
    otherwise
        error('Filename must be specified as a string or character argument.')
end

%% Read file contents
idx = 0;
while true
    str = fgetl(fID);
    if str == -1
        break
    end

    idx = idx+1;
    mat = regexp(str, '(\w*)', 'match');
    
    if isempty(mat)
        fprintf('No Variable Name\n\tLINE %d\n\t%s\n',idx,str);
        continue
    end

    eval( sprintf('vars.%s',str) );
end

%% Close file 
fclose(fID);
