function [ZERO,fast,cellOut] = parseVarargin_ZERO_fast(cellIn)
% PARSEVARARGIN_ZERO_FAST parses a variable input argument for ZERO and
% fast parameters.
%   [ZERO,fast] = parseVarargin_ZERO_fast(cellIn)
%   
%   % TODO - Update documentation and finish function.
%
%   Output(s)
%       ZERO - [OPTIONAL] positive value that is sufficiently close to zero
%              or assumed zero (e.g. ZERO = 1e-8). If ZERO is not
%              specified, a default value is used.
%       fast - [OPTIONAL] true/false logical value indicating whether to
%              skip checking SE(N). Choosing fast = true ignores specified
%              ZERO.
%                fast = true    - Skip checking if H \in SE(N)
%                fast = [false] - Check if H \in SE(N)