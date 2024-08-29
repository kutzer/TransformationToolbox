function strTex = str2LaTeX(str)
% STR2TEX converts a character array into a formatted LaTeX equation
%
%   strTex = str2LaTeX(str)
%
%   Input(s)
%       str - character array containing the equation to be converted
%
%   Output(s)
%       strTex - character array formatted for LaTeX
%
%   NOTES:
%       (1) Superscripts and subscripts can be added using H_aTob
%           convention. For example:
%               H_1To0 --> H_{1}^{0}
%               H_To0  --> H^{0}
%               R_aToStar --> R_{a}^{*}
%       (2) Modifiers such as bar, underline, vec, etc can be added using
%           camel case text prior to the variable. For example:
%               vecUnderlineB_To6 --> \vec{\underline{B}}^{6}
%
%   See also sym2LaTeX array2LaTeX
%
%   M. Kutzer, 07Feb2024, USNA

% TODO - support partial derivatives
% TODO - allow eta to be replaced
% etc.

%% Check input(s)
% TODO - check inputs 

%% Convert string
% Initialize strTex
strTex = str;

% Replace known LaTeX math terms
% -> Define term(s)
symTerm = {...
    '*',... % ----------------> * must come before Star
    'sin','cos','tan','asin','acos','atan','log',...
    'logm','Star',...
    'alpha','beta','gamma','delta','epsilon','zeta','theta',...
    'iota','kappa','lambda','mu','nu','xi','pi','rho','sigma',...
    'tau','upsilon','chi','psi','omega','phi'};%,...
%'eta'}; % ----------------> eta is not currently supported.
% -> Define replacement(s)
texTerm = {...
    ' ',... % ----------------> * must come before Star
    '\sin','\cos','\tan','\asin','\acos','\atan','\log',...
    '\log','*',...
    '\alpha','\beta','\gamma','\delta','\epsilon','\zeta','\theta',...
    '\iota','\kappa','\lambda','\mu','\nu','\xi','\pi','\rho','\sigma',...
    '\tau','\upsilon','\chi','\psi','\omega','\phi'};%,...
%'\eta'}; % ----------------> eta is not currently supported.
% -> Perform the replacements
for i = 1:numel(texTerm)
    strTex = strrep(strTex, symTerm{i}, texTerm{i});
end

% Replace superscript & subscript terms
% -> Define the pattern
pattern = '([A-Za-z\d]+)_([A-Za-z\d\*]+)To([A-Za-z\d\*]+)';
% -> Define replacement
replacement = '$1_{$2}^{$3}';
% -> Perform the replacement using regular expression
strTex = regexprep(strTex, pattern, replacement);

% Replace subscript only terms
% -> Define the pattern
pattern = '([A-Za-z\d]+)_([A-Za-z\d]+)';
% -> Define replacement
replacement = '$1_{$2}';
% -> Perform the replacement using regular expression
strTex = regexprep(strTex, pattern, replacement);

% Replace superscript only terms
% -> Define the pattern
pattern = '([A-Za-z\d]+)_To([A-Za-z\d\*]+)';
% -> Define replacement
replacement = '$1^{$2}';
% -> Perform the replacement using regular expression
strTex = regexprep(strTex, pattern, replacement);

% Replace specific modifier terms
symTerm = {'vec','underline','bar','hat'};
for i = 1:numel(symTerm)
    % -> Define the pattern
    pattern = sprintf('(%s)([A-Za-z\\d]+)',symTerm{i});
    % -> Define replacement
    replacement = '\\$1{$2}';
    % -> Perform the replacement using regular expression
    strTex = regexprep(strTex, pattern, replacement);
end

% Account for special cases (sqrt)
% TODO - This doesn't work!
%{
% -> Define the pattern
pattern = '\((.*?)\)\^\(1/2\)';
% -> Define the replacement
replacement = '\\sqrt{$1}';
% -> Perform the replacement using regular expression
strTex = regexprep(strTex, pattern, replacement);
%}

% Account for fractions (this should not be the case)
% TODO - This doesn't work!
%{
% -> Define the pattern
pattern = '([^/]+)/([^/]+)';
% -> Define the replacement
replacement = '\\frac{$1}{$2}';
% -> Perform the replacement using regular expression
strTex = regexprep(strTex, pattern, replacement);
%}

% Replace remaining terms
% -> Define term(s)
symTerm = {'(',')','[',']'};
% -> Define replacement(s)
texTerm = {'\left(','\right)','\left[','\right]'};
% -> Perform the replacements
for i = 1:numel(texTerm)
    strTex = strrep(strTex, symTerm{i}, texTerm{i});
end
