function kStr = transformsStr2kinematicStr(tforms)
% TRANSFORMSSTR2KINEMATICSTR compiles individual transformations into a
% string of transformations
%   kStr = TRANSFORMSSTR2KINEMATICSTR(tforms)
%
%   Input(s)
%       tforms - cell array containing an ordered set of individual motion
%                primitives (Tx, Ty, Tz, Rx, Ry, and Rz) contained in the 
%                kinematic string. 
%
%   Output(s)
%       str - string of kinematic motion primitives (Tx, Ty, Tz, Rx, Ry,
%             and Rz) each containing fixed values or variables. 
%
%   M. Kutzer, 06Jul2021, USNA

%% Check inputs
narginchk(1,1);
% TODO - check input string

%% Compile string
kStr = tforms{1};
for i = 2:numel(tforms)
    kStr = sprintf('%s*%s',kStr,tforms{i});
end