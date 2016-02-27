function TransformationToolboxUpdate
% TRANSFORMATIONTOOLBOXUPDATE download and update the Transformation 
% Toolbox. 
%
%   M. Kutzer 27Feb2016, USNA

% TODO - Find a location for "TransformationToolbox Example SCRIPTS"
% TODO - update function for general operation

%% Check current version
A = TransformationToolboxVer;

%% Setup temporary file directory
fprintf('Downloading the Transformation Toolbox...');
tmpFolder = 'TransformationToolbox';
pname = fullfile(tempdir,tmpFolder);

%% Download and unzip toolbox (GitHub)
url = 'https://github.com/kutzer/TransformationToolbox/archive/master.zip';
try
    fnames = unzip(url,pname);
    fprintf('SUCCESS\n');
    confirm = true;
catch
    confirm = false;
end

%% Check for successful download
if ~confirm
    error('Failed to download updated version of Transformation Toolbox.');
end

%% Find base directory
install_pos = strfind(fnames,'installTransformationToolbox.m');
sIdx = cell2mat( install_pos );
cIdx = ~cell2mat( cellfun(@isempty,install_pos,'UniformOutput',0) );

pname_star = fnames{cIdx}(1:sIdx-1);

%% Get current directory and temporarily change path
cpath = cd;
cd(pname_star);

%% Install ScorBot Toolbox
installTransformationToolbox(true);

%% Move back to current directory and remove temp file
cd(cpath);
[ok,msg] = rmdir(pname,'s');
if ~ok
    warning('Unable to remove temporary download folder. %s',msg);
end

%% Complete installation
fprintf('Installation complete.\n');

end
