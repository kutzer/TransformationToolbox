function [bin,msg] = isSkewSymmetric(M,ZERO)
% ISSKEWSYMMETRIC checks a matrix to see if it is skew-symmetric
%   bin = ISSKEWSYMMETRIC(M) checks nxn matrix "M" for skew-symmetry. If 
%   the M is skew-symmetric, this function returns "1", "0" is returned 
%   otherwise.
%
%   bin = ISSKEWSYMMETRIC(M,ZERO) allows a user to specify an acceptable
%   zero value.
%
%   [bin,msg] = ISSKEWSYMMETRIC(___) returns a description of why the 
%   matrix is not skew symmetric.
%
%   M. Kutzer 10Oct2014, USNA

% Updates
%   24Jan2017 - Updated to return a binary and associated message
%   16Jan2018 - Updated to define ZERO based on M.
%   07Feb2018 - Updated to refine ZERO estimate based on specific test
%               condition.
%   18Nov2021 - Added optional ZERO input

%% Check input(s)
narginchk(1,2);
if nargin < 2
    ZERO = [];
end

%% Define ZERO scaling term
ZERO_scale = 1e1;

%% Check dimensions of M
if size(M,1) ~= size(M,2) || ~ismatrix(M)
    bin = false;
    msg = '"M" is not a square matrix.';
    return
end

%% Check for custom functions
%TODO - check for zeroFPError.m

%% Check for skew-symmetric
msg = [];
bin = true; % assume skew-symmetric matrix

if strcmpi(class(M),'double') || strcmpi(class(M),'single')
    %ZERO = 50*eps( max(abs(reshape(M,1,[]))) );
    if ~isempty(ZERO)
        ZERO_i = ZERO_scale * max( reshape(eps(M),1,[]) );
    else
        ZERO_i = ZERO;
    end
    chk = zeroFPError( abs(M + transpose(M)),ZERO_i );
else
    chk = zeroFPError( abs(M + transpose(M)) );
end

try
    % Check if term contains any symbolic variables
    logical(chk);
catch
    % Try simplifying complicated term one more time
    chk = zeroFPError(chk,ZERO);
end

try 
    if max( reshape(chk,1,[]) ) > 0  % converts "chk" to logical, it will 
                                     % throw an error if for symbolic 
                                     % arguments are still in the 
                                     % expression.
        bin = false;
        msg = sprintf(['"abs(M + transpose(M))" returned at least one\n\t\t',...
            'non-zero value: %.15f.'], max( reshape(chk,1,[]) ));
        return
    end
catch
    bin = false;
    msg = 'Unable to check for skew-symmetry.';
    return
end

