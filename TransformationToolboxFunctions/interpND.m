function v_out = interpND(v0,v1,delta)
% INTERPND creates a linear interpolation between two vectors given a
% specified delta.
%   v_out = interpND(v0,v1,delta)
%
%   Input(s)
%       v0    - Nx1 array specifying the initial N-D configuration
%       v1    - Nx1 array specifying the final N-D configuration 
%       delta - [OPTIONAL] positive scalar value specifying the maximum 
%               change in any given value. If unspecified, delta = 0.01.
%
%   Output(s)
%       v_out - NxM array specifying the linear interpolation from v0 to v1
%               with a maximum step size of approximately delta for any 
%               given dimension. Note that the actual minimum change will
%               be less than or equal to delta depending on the minimum
%               spacing between v0 and v1.
%
%   M. Kutzer, 13Oct2022, USNA

%% Check input(s)
narginchk(2,3);

if nargin < 3
    delta = 0.01;
end

if numel(v0) ~= numel(v1)
    error('v0 and v1 must be the same dimension.');
end

if numel(delta) ~= 1
    error('delta must be specified as a positive scalar value.');
end

if delta <= 0
    error('delta must be specified as a positive scalar value.');
end

N = numel(v0);
if numel(v1) ~= N
    warning('Inputs v0 and v1 should be %dx1 arrays, output will be %dxM.',N,N);
end

%% Find the minimum spacing
v0 = reshape(v0,[],1);
v1 = reshape(v1,[],1);

dv = max( abs(v1 - v0) );

interp_dv = 0:delta:dv;
M = numel(interp_dv);

if interp_dv(end) < dv
    M = M+1;
end

%% Interpolate
v_out = zeros(N,M);
for k = 1:N
    v_out(k,:) = linspace(v0(k),v1(k),M);
end