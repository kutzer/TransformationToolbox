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

A.Name = 'Transformation Toolbox';
A.Version = '1.0.8';
A.Release = '(R2019a)';
A.Date = '17-Oct-2019';
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