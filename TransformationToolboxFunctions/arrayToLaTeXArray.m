function str = arrayToLaTeXArray(M,dec,lrc)
% ARRAYTOLATEXARRAY creates a latex array string given a matrix and
% specificied number decimal places
%   str = arrayToLaTeXArray(M)
%
%   str = arrayToLaTeXArray(M,dec)
%
%    str = arrayToLaTeXArray(M,dec,alignment)
%
%   Input(s)
%         M - nxm array containing numeric and/or symbolic values
%       dec - [OPTIONAL] specifies the number of decimal places for
%             specified values. The default value for dec is 4.
%       lrc - [OPTIONAL] specifies the alignment within the array
%             using an element of {'l',['c'],'r'}.
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
%   See also array2LaTeX sym2LaTeX str2LaTeX
%
%   M. Kutzer, 27Oct2021, USNA

% Updates
%   11May2022 - Function completed
%   05Feb2024 - Added check dec and alignment
%   05Feb2024 - Added symbolic support
%% Check input(s)
narginchk(1,3);

if nargin < 3
    lrc = 'c';
end

if nargin < 2
    dec = 4;
end

% TODO - check M

% Check dec
if dec < 0 || numel(dec) ~= 1 || round(dec) ~= dec
    error('Number of decimals must be a scalar integer that is greater than or equal to 0.');
end

% Check alignment (use contains)
lrc = lower(lrc);
supportedAlignment = {'l','c','r'};
if ~any( matches(supportedAlignment,lrc) )
    error('Alignment must be ''l'', ''c'', or ''r''');
end

%% Define decimal string
if ~strcmpi( class(M), 'sym')
    maxVal = max(reshape(M,1,[]));
    nValsLHS = numel( num2str(round(maxVal)) );
    nValsRHS = dec;
    nVals = nValsLHS + nValsRHS + 2; % account for decimal and negative sign
    decStr = ['%',num2str(nVals),'.',num2str(dec),'f'];
end

%% Display array
str = sprintf('\\left( \\begin{array}{%s}\n',repmat(lrc,1,size(M,2)));
for i = 1:size(M,1)
    if i > 1
        str = sprintf('%s \\\\\n',str);
    end
    for j = 1:size(M,2)
        if j < size(M,2)
            if ~strcmpi( class(M), 'sym')
                % numeric inputs
                str = sprintf(['%s',decStr,' & '],str,M(i,j));
            else
                % apply for symbolic inputs
                str = sprintf('%s%s & ',str,sym2LaTeX(M(i,j),dec));
            end
        else
            if ~strcmpi( class(M), 'sym')
                % numeric inputs
                str = sprintf(['%s',decStr],str,M(i,j));
            else
                % apply for symbolic inputs
                str = sprintf('%s%s',str,sym2LaTeX(M(i,j),dec));
            end
        end
    end
end
str = sprintf('%s \\\\\n\\end{array}\\right)\n',str);

fprintf('%s',str);