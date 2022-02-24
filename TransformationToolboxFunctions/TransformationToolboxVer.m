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

A.Name = 'Transformation Toolbox';
A.Version = '1.1.10';
A.Release = '(R2019b)';
A.Date = '24-Feb-2022';
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