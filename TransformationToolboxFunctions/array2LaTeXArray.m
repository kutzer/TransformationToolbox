function str = array2LaTeXArray(varargin)
% ARRAY2LATEXARRAY creates a latex array string given a matrix and
% specificied number decimal places
%   str = array2LaTeXArray(M)
%
%   str = array2LaTeXArray(M,dec)
%
%    str = array2LaTeXArray(M,dec,lrc)
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
%   See also sym2LaTeX str2LaTeX
%
%   M. Kutzer, 07Feb2024, USNA

% NOTE - This function is renaming arrayToLaTeXArray.m
str = arrayToLaTeXArray(varargin{:});