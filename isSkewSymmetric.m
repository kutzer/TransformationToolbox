function bin= isSkewSymmetric(M)
% isSkewSymmetric checks a matrix to see if it is skew-symmetric
%   isSkewSymmetric(M) checks nxn matrix "M" for skew-symmetry. If the M is
%   skew-symmetric, this function returns "1", "0" is returned otherwise.
%
%   (c) M. Kutzer 10Oct2014, USNA

%% Check dimensions of M
if size(M,1) ~= size(M,2) || ~ismatrix(M)
    error('"M" must be a square matrix.');
end

%% Check for custom functions
%TODO - check for zeroFPError.m

%% Check for skew-symmetric
bin = true; % assume skew-symmetric matrix

chk = zeroFPError( sum(sum( abs(M + transpose(M)) )) );
try
    % Check if term contains any symbolic variables
    logical(chk);
catch
    % Try simplifying complicated term one more time
    chk = zeroFPError(chk);
end

try 
    if abs(chk) > 0  % converts "chk" to logical, it will throw an error if  
                     % for symbolic arguments are still in the expression.
        bin = false;
        return
    end
catch
    bin = false;
    return
end

