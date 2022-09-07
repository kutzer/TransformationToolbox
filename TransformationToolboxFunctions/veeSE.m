function v = veeSE(h,varargin)
% VEESE converts an NxN element of se(n) (the Lie algebra associated with
% SE(n), sometimes referred to as "little se") into a vector defined by the
% basis elements of se(N).
%   v = VEESE(M) converts an NxN matrix "h" into an Mx1 vector "v".
%
%   r = VEESE(___,ZERO)
%
%   r = VEESE(___,fast)
%
%   v = VEESE(___,'fast') <--- !!! LEGACY SYNTAX !!!
%
%   Input(s)
%       r - NxN element of se(N)
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not   
%              specified, a default value is used.
%       fast - [OPTIONAL] true/false logical value indicating whether to
%              skip checking for skew-symmetry. Choosing fast = true 
%              ignores specified ZERO. 
%                fast = true    - Skip checking if h is skew symmetric
%                fast = [false] - Check if h is skew symmetric
%
%   Output(s)
%       v - Mx1 vector containing the unique values of h ordered using 
%       3 x 3 matrix -> $v \in \mathbb{R}^3$ (2D rigid body motion)
%       4 x 4 matrix -> $v \in \mathbb{R}^6$ (3D rigid body motion)
%       5 x 5 matrix -> $v \in \mathbb{R}^10$
%       6 x 6 matrix -> $v \in \mathbb{R}^15$
%       ...
%       N x N matrix -> $v \in \mathbb{R}^M$
%
%   See also wedgeSE, seBasis, veeSO, wedgeSO.
%
%   M. Kutzer, 04Jan2017, USNA

% Updates
%   26Jan2022 - Added increased "zero" value of 1e-8
%   06Sep2022 - Updated to include ZERO and "standardized" fast optional 
%               inputs

%% Default options
ZERO = 1e-8;
fast = false;

%% Default options
narginchk(1,3);

if nargin > 1
    for i = 1:numel(varargin)
        switch lower( class(varargin{i} ))
            case 'char'
                % Legacy input
                str = varargin{i};
                switch lower(str)
                    case 'fast'
                        fast = true;
                end
            case 'string'
                % Legacy input
                str = char( varargin{i} );
                switch lower(str)
                    case 'fast'
                        fast = true;
                end
            case 'logical'
                fast = varargin{i};
            otherwise
                if numel(varargin{i}) ~= 1
                    error('Numeric optional inputs must be scalar values.');
                end

                if varargin{i} == 0 || varargin{i} == 1
                    fast = logical(varargin{i});
                else
                    ZERO = varargin{i};
                end
        end
    end
end

% Check zero
if ZERO < 0
    error('ZERO value must be greater or equal to 0.')
end

%% Check M
N = size(h,1);
if ~fast
    if exist('isSkewSymmetric','file')
        [bin,msg] = isSkewSymmetric(h(1:(N-1),1:(N-1)),ZERO);
        if ~bin
            warning('"M(1:(N-1),1:(N-1))" must be skew-symmetric.\n\t -> %s',msg);
        end
    else
        warning('The function "isSkewSymmetric" was not found. Skipping skew-symmetry check.');
    end
end

%% $M \in se(N)$
n = N-1;
e = seBasis(n);
m = numel(e);
for idx = 1:m
    [i,j] = find(e{idx} == 1);
    v(idx,1) = h(i,j);
end