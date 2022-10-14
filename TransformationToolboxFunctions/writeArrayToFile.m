function varargout = writeArrayToFile(fname,var,i)
% WRITEARRAYTOFILE writes an array to a specified filename or file ID.
%   writeArrayToFile(fname,var,i)
%
%   writeArrayToFile(fID,var,i)
%
%   fID = writeArrayToFile(fname,var,i)
%
%   Input(s)
%       fname - character array specifying filename for writing 
%       fID   - file ID (e.g. "fID = fopen(fname,'a')") 
%       var   - array to be written to file
%       i     - [OPTIONAL] positive integer specifying variable index
%
%   Output(s)
%       fID - file ID of file being written. 
% 
%   NOTE(s):
%       The file is left open if:
%           (1) a file ID is provided as an input
%           (2) a filename is provided as an input, and no output is
%               defined      
%
%   Example (add a single variable)
%       fname = 'EXAMPLE_writeArrayToFile.txt';
%       X_tst = rand(1,2,2);
%       writeArrayToFile(fname,X_tst);
%
%       'EXAMPLE_writeArrayToFile.txt' File Contents:
%           X_tst = reshape([___, ___, ___, ___],1,2,2);\n
%
%   Example (add multiple variables)
%       fname = 'EXAMPLE_writeArrayToFile.txt';
%       fID = fopen(fname,'a');
%       for i = 1:3
%           H_a2b = rand(2,2);
%           writeArrayToFile(fID,H_a2b,i);
%           H_b2c = rand(2,2);
%           writeArrayToFile(fID,H_b2c,i);
%       end
%       fclose(fID);
%
%       'EXAMPLE_writeArrayToFile.txt' File Contents:
%           H_a2b{1} = reshape([___, ___, ___, ___],2,2);\n
%           H_b2c{1} = reshape([___, ___, ___, ___],2,2);\n
%           H_a2b{2} = reshape([___, ___, ___, ___],2,2);\n
%           H_b2c{2} = reshape([___, ___, ___, ___],2,2);\n
%           H_a2b{3} = reshape([___, ___, ___, ___],2,2);\n
%           H_b2c{3} = reshape([___, ___, ___, ___],2,2);\n
%
%   See also readFileToArray
%
%   M. Kutzer, 14Oct2022, USNA

%% Check input(s)
narginchk(2,3);

tf_fclose = false;
switch class(fname)
    case {'char','string'}
        % Open file
        fID = fopen(fname,'a');
        if fID < 0 
            error('Unable to open "%s". Try a new filename.',fname);
        end

        % Update close file flag
        if nargout == 0
            tf_fclose = true;
        end

    otherwise
        %error('Filename must be specified as a string or character argument.')
        fID = fname;

end

if nargin < 3
    i = [];
end
%% Get input variable name
varName = inputname(2);

%% Create string for file
% Create initial string
if ~isempty(i)
    str_0 = sprintf('%s{%d} = reshape([',varName,i);
else
    str_0 = sprintf('%s = reshape([',varName);
end

% Create numeric values in string
varRSHP = reshape(var,1,[]);
if numel(varRSHP) > 1
    str_i = sprintf('%.6f,',varRSHP(1:end-1));
end
str_i = sprintf('%s%.6f],',str_i,varRSHP(end));

% Define reshape values 
n = ndims(var);
N = zeros(1,n);
for j = 1:n
    N(j) = size(var,j);
end

if numel(N) > 1
    str_j = sprintf('%d,',N(1:end-1));
end
str_j = sprintf('%s%d);',str_j,N(end));

% Combine string
str = sprintf('%s%s%s',str_0,str_i,str_j);

%% Show string to be written to file
%fprintf('%s\n',str);

%% Write to file
fprintf(fID,'%s\n',str);

%% Package output(s)
if nargout > 0
    varargout{1} = fID;
else
    varargout = {};
end

%% Close file (if applicable)
if tf_fclose
    fclose(fID);
end