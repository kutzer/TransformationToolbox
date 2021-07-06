function tforms = kinematicStr2transformsStr(kStr)
% KINEMATICSTR2TRANSFORMSSTR isolates individual transforms from a string of
% transformations
%   tforms = KINEMATICSTR2TRANSFORMSSTR(kStr)
%
%   Input(s)
%       str - string of kinematic motion primitives (Tx, Ty, Tz, Rx, Ry,
%             and Rz) each containing fixed values or variables. 
%
%   Output(s)
%       tforms - cell array containing an ordered set of individual motion
%                primitives (Tx, Ty, Tz, Rx, Ry, and Rz) contained in the 
%                kinematic string. 
%
%   M. Kutzer, 06Jul2021, USNA

%% Check inputs
narginchk(1,1);
% TODO - check input string
 
%% Parse string
idx = findstr(kStr,'*');
idx = [0,idx,numel(kStr)+1]; % Append "first" and "last" index
for i = 2:numel(idx)
    tforms{i-1} = kStr( (idx(i-1)+1):(idx(i)-1) );
end
