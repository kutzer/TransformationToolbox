function str = arrayToLaTeXArray(M,dec)
% ARRAYTOLATEXARRAY creates a latex array string given a matrix and
% specificied number decimal places
%   str = arrayToLaTeXArray(M)
%
%   str = arrayToLaTeXArray(M,dec)
%
%   Input(s)
%       M   - nxm array containing numeric and/or symbolic values
%       dec - [OPTIONAL] specifies the number of decimal places for
%             specified values. The default value for dec is 4.
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

error('This function is incomplete')