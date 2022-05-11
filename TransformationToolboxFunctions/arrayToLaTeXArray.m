function str = arrayToLaTeXArray(M,dec,alignment)
% ARRAYTOLATEXARRAY creates a latex array string given a matrix and
% specificied number decimal places
%   str = arrayToLaTeXArray(M)
%
%   str = arrayToLaTeXArray(M,dec)
%
%    str = arrayToLaTeXArray(M,dec,alignment)
%
%   Input(s)
%               M - nxm array containing numeric and/or symbolic values
%             dec - [OPTIONAL] specifies the number of decimal places for
%                   specified values. The default value for dec is 4.
%       alignment - [OPTIONAL] specifies the alignment within the array
%                   using an element of {'l',['c'],'r'}.
%
%   Output(s)
%       str - string containing LaTeX array format
%
%   Example
%       M = [1.0052, 0; -0, 5002.3]
%
%       str = arrayToLaTeXArray(M,5)
%
%       str = '...
%           1.00520 & 0       \\ \n
%           0       & 5.002.3'
%
%   M. Kutzer, 27Oct2021, USNA

% Updates
%   11May2022 - Function completed

%% Check input(s)
narginchk(1,3);

if nargin < 3
    alignment = 'c';
end

if nargin < 2
    dec = 4;
end

% TODO - check M
% TODO - check dec
% TODO - check alignment (use contains)

%% Define decimal string
maxVal = max(reshape(M,1,[]));
nValsLHS = numel( num2str(round(maxVal)) );
nValsRHS = dec;
nVals = nValsLHS + nValsRHS + 2; % account for decimal and negative sign
decStr = ['%',num2str(nVals),'.',num2str(dec),'f'];

%% Display array
str = sprintf('\\left( \\begin{array}{%s}\n',repmat(alignment,1,size(M,2)));
for i = 1:size(M,1)
    if i > 1
        str = sprintf('%s \\\\\n',str);
    end
    for j = 1:size(M,2)
        if j < size(M,2)
            str = sprintf(['%s',decStr,' & '],str,M(i,j));
        else
            str = sprintf(['%s',decStr],str,M(i,j));
        end
    end
end
str = sprintf('%s \\\\\n\\end{array}\\right)\n',str);

fprintf('%s',str);