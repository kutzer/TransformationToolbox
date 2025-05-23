function varargout = TransformationToolboxVer
% TRANSFORMATIONTOOLBOXVER displays the Transformation Toolbox information.
%   TRANSFORMATIONTOOLBOXVER displays the information to the command prompt.
%
%   A = TRANSFORMATIONTOOLBOXVER returns in A the sorted struct array of  
%   version information for the Transformation Toolbox.
%     The definition of struct A is:
%             A.Name      : toolbox name
%             A.Version   : toolbox version number
%             A.Release   : toolbox release string
%             A.Date      : toolbox release date
%
%   M. Kutzer 27Feb2016, USNA

% Updates
%   07Mar2018 - Updated to include try/catch for required toolbox
%               installations
%   15Mar2018 - Updated to include msgbox warning when download fails
%   17Oct2019 - Added stereoCorrespondenceToPoint.m
%   13Jul2020 - Updated setTriad 'scale' parameter
%   25Sep2020 - Updated to add nearestSO and nearestSE
%   08Jan2021 - Updated ToolboxUpdate
%   14Apr2021 - Added solveAXeqXBinSE.m
%   19Apr2021 - Updated covSE and meanSE to use invSE
%   29Jun2021 - Updated to add rigidBodyTree support
%   02Nov2021 - Updated meanSE with forced real skew-symmetry 
%   18Nov2021 - Updated meanSE, isSO, isSE, and isZERO to include a
%               user-specified optional ZERO input
%   26Jan2022 - Updated SO/SE log/exp functions
%   26Jan2022 - Updated default ZERO values
%   27Jan2022 - Updated to include SEtoSO
%   24Feb2022 - Updated calculateJacobian
%   24Feb2022 - Updated nCross 
%   28Feb2022 - Added comments to saved function from calculateJacobian
%   03Mar2022 - Added limitVector and numericIkin functions
%   03Mar2022 - Added showStatus option to calculateJacobian
%   25Apr2022 - Added AxesLabels "property" to setTriad
%   25Aug2022 - Accounted for NaN values in isZero
%   02Sep2022 - Removed SE(3) assumption and added coupled/decoupled
%               options to meanSE
%   07Sep2022 - Updated covSE, veeSE, invSE, etc to create common approach
%   27Sep2022 - Added multiviewPointsToWorld
%   14Oct2022 - Added writeArrayToFile and readArrayToFile
%   25Oct2022 - Preallocate memory and show waitbar for large data sets in 
%               solveAXeqXBinSE
%   23Aug2024 - Improved triad support
%   31Jan2025 - Added revolute and prismatic joint visualization
%   02Apr2025 - Updated to make "bin" always logical in isSO and isSE
%   22May2025 - Updated for local user install

A.Name = 'Transformation Toolbox';
A.Version = '1.1.23';
A.Release = '(R2019b)';
A.Date = '23-May-2025';
A.URLVer = 1;

msg{1} = sprintf('MATLAB %s Version: %s %s',A.Name, A.Version, A.Release);
msg{2} = sprintf('Release Date: %s',A.Date);

n = 0;
for i = 1:numel(msg)
    n = max( [n,numel(msg{i})] );
end

fprintf('%s\n',repmat('-',1,n));
for i = 1:numel(msg)
    fprintf('%s\n',msg{i});
end
fprintf('%s\n',repmat('-',1,n));

if nargout == 1
    varargout{1} = A;
end