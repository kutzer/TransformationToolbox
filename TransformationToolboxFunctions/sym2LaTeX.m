function str = sym2LaTeX(val,dec)
% SYM2LATEX converts a symbolic scalar term into a latex equation string.
%
%   str = sym2LaTeX(val,dec)
%
%   Input(s)
%       val - scalar symbolic expression
%       dec - [OPTIONAL] number of decimals to include in numeric values
%
%   Output(s)
%       str - string representing scalar symbolic expression
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
%   See also array2LaTeXArray str2LaTeX
%
%   M. Kutzer, 05Feb2024, USNA

%% Check input(s)

% TODO - check "val"

if nargin < 2
    dec = 4;
end

% Check dec
if dec < 0 || numel(dec) ~= 1 || round(dec) ~= dec
    error('Number of decimals must be a scalar integer that is greater than or equal to 0.');
end

%% Define string for printing decimals
decStr = ['%','.',num2str(dec),'f'];

%% Convert terms
str = '';
if strcmpi( class(val), 'sym') % apply for symbolic inputs

    % Isolate numerator and denominator
    [n,d] = numden(val); % account for expressions that are fractions
    
    if n == 0
        % Account for special case (numerator of 0)
        str = sprintf(decStr,n);
    else
        % Account for standard case 

        % Numerator
        [c,s] = coeffs(n);
        if ~isempty(c) && ~isempty(s)
            for k = 1:numel(c)
                % Convert symbolic expression to latex
                [~,strPos,strSgn] = cs2str(c(k),s(k),decStr);

                % Append string
                str = appendNewString(str,strPos,strSgn);
            end
        else
            n = 0;
        end

        % Denominator
        if d ~= 1
            strD = '';
            [c,s] = coeffs(d);
            for k = 1:numel(c)
                % Convert symbolic expression to latex
                [~,strPos,strSgn] = cs2str(c(k),s(k),decStr);

                % Append string
                strD = appendNewString(strD,strPos,strSgn);
            end
            str = sprintf('\\frac{%s}{%s}',str,strD);
        end

    end
else

    str = sprintf(decStr,val);

end

%% Account for values of 0 and 1
strZero = sprintf(decStr,0);
if matches(str,strZero)
    str = '0';
end

strOne  = sprintf(decStr,1);
if matches(str,strOne)
    str = '1';
end

end

%% ------------------------- Internal Functions ---------------------------

% -------------------------------------------------------------------------
function [strNew,strPos,strSgn] = cs2str(ci,si,decStr)


% Define coeffient
if abs(ci) ~= 1
    strNew = sprintf(decStr,ci);
    strPos = sprintf(decStr,abs(ci));
else
    strNew = '';
    strPos = '';
end
strSgn = sign(ci);

% Convert symbolic term to string
strSym = char(si);

% Account for special case (si == 1)
if matches(strSym,'1')
    strNew = sprintf(decStr,ci);
    strPos = sprintf(decStr,abs(ci));
    return
end

% Process symbolic string
strSym = str2LaTeX(strSym);

% Combine strings
strNew = sprintf('%s%s',strNew,strSym);
strPos = sprintf('%s%s',strPos,strSym);

% Account for empty strSym
if isempty(strSym)
    strNew = sprintf(decStr,ci);
    strPos = sprintf(decStr,abs(ci));
end

end

% -------------------------------------------------------------------------
function str = appendNewString(str,strPos,strSgn)

if isempty(str)
    if strSgn < 0
        str = sprintf('-%s',strPos);
    else
        str = strPos;
    end
else
    if strSgn < 0
        str = sprintf('%s - %s',str,strPos);
    else
        str = sprintf('%s + %s',str,strPos);
    end
end

end
