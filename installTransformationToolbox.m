function installTransformationToolbox(replaceExisting,skipAdmin)
% INSTALLTRANSFORMATIONTOOLBOX installs Transformation Toolbox for MATLAB.
%   INSTALLTRANSFORMATIONTOOLBOX installs Transformation Toolbox into the following 
%   locations:
%                        Source: Destination
%     TransformationToolboxFunctions: matlabroot\toolbox\transformation
%       TransformationToolboxSupport: matlabroot\toolbox\transformation\TransformationToolboxSupport 
%
%   INSTALLTRANSFORMATIONTOOLBOX(true) installs Transformation Toolbox regardless of
%   whether a copy of the Transformation toolbox exists in the MATLAB root.
%
%   INSTALLTRANSFORMATIONTOOLBOX(false) installs Transformation Toolbox only if no copy 
%   of the Transformation toolbox exists in the MATLAB root.
%
%   M. Kutzer 17Feb2016, USNA

% Updates
%   08Jan2021 - Corrected questdlg
%   22May2025 - Enable local user installation

%% Assign tool/toolbox specific parameters
dirName = 'transformation';
toolboxContent  = 'TransformationToolboxFunctions';
toolboxExamples = 'TransformationToolbox Example SCRIPTS';
toolboxName = 'Transformation Toolbox';
toolboxShort = strrep(toolboxName, ' ', '');

%% Define toolbox directory options
toolboxPathAdmin = fullfile(matlabroot,'toolbox',dirName);
toolboxPathLocal = fullfile(prefdir,'toolbox',dirName);
toolboxPathExamples = fullfile(userpath,sprintf('Examples, %s',toolboxName));

%% Check if folders exist
isPathAdmin = isfolder(toolboxPathAdmin);
isPathLocal = isfolder(toolboxPathLocal);
isPathExamples = isfolder(toolboxPathExamples);

%% Check if folders are in the MATLAB path
allPaths = path;
allPaths = strsplit(allPaths,pathsep);

inPathAdmin = any(matches(allPaths,toolboxPathAdmin),'all');
inPathLocal = any(matches(allPaths,toolboxPathLocal),'all');
inPathExamples = any(matches(allPaths,toolboxPathExamples),'all');

%% Check inputs
if nargin < 1
    replaceExisting = [];
end
if nargin < 2
    skipAdmin = false;
end

%% Check for admin write access
isAdmin = checkWriteAccess(matlabroot);
isLocal = checkWriteAccess(prefdir);
isExample = checkWriteAccess(userpath);

%% Check for basic write access
if ~isLocal
    warning('No local write access?');
    return
end

%% Installation error solution(s)
adminSolution = sprintf(...
    ['Possible solution:\n',...
     '\t(1) Close current instance of MATLAB\n',...
     '\t(2) Open a new instance of MATLAB "as administrator"\n',...
     '\t\t(a) Locate MATLAB shortcut\n',...
     '\t\t(b) Right click\n',...
     '\t\t(c) Select "Run as administrator"\n']);

%% Prompt user if not an admin
if ~isAdmin && ~skipAdmin
    choice = questdlg(sprintf(...
        ['MATLAB is running without administrative privileges.\n',...
        'Would you like to install the %s locally?'],toolboxName),...
        sprintf('Local Install %s',toolboxName),...
        'Yes','No','Cancel','Yes');
    % Replace existing or cancel installation
    switch choice
        case 'Yes'
            skipAdmin = true;
        case 'No'
            fprintf('Unable to write perform installation.\n\n');
            fprintf('To install as an administrator - %s',adminSolution);
            return
        case 'Cancel'
            fprintf('Action cancelled.\n');
            return
        otherwise
            error('Unexpected response.');
    end
end

%% Check for existing toolbox
if skipAdmin
    isToolbox = isPathLocal;
    toolboxPath = toolboxPathLocal;
else
    isToolbox = isPathAdmin;
    toolboxPath = toolboxPathAdmin;
end

%% Check for toolbox directory
if isToolbox
    % Apply replaceExisting argument
    if isempty(replaceExisting)
        choice = questdlg(sprintf(...
            ['MATLAB path already contains the %s.\n',...
            'Would you like to replace the existing toolbox?'],toolboxName),...
            sprintf('Replace Existing %s',toolboxName),...
            'Yes','No','Cancel','Yes');
    elseif replaceExisting
        choice = 'Yes';
    else
        choice = 'No';
    end
    % Replace existing or cancel installation
    switch choice
        case 'Yes'
            % Remove existing paths
            removePath(toolboxName,...
                toolboxPathAdmin,inPathAdmin,isPathAdmin,isAdmin);
            removePath(toolboxName,...
                toolboxPathLocal,inPathLocal,isPathLocal,isLocal);
            removePath(toolboxName,...
                toolboxPathExamples,inPathExamples,inPathExamples,isExample);
        case 'No'
            fprintf('%s currently exists, installation cancelled.\n',toolboxName);
            return
        case 'Cancel'
            fprintf('Action cancelled.\n');
            return
        otherwise
            error('Unexpected response.');
    end
end

%% Migrate toolbox folder contents
% Toolbox contents
migrateContent(toolboxContent,toolboxPath,toolboxName);
% Example files
try
migrateContent(toolboxExamples,toolboxPathExamples,...
    sprintf('%s Examples',toolboxName));
catch
    fprintf('Unable to migrate examples.');
end

%% Save toolbox path
%addpath(genpath(toolboxRoot),'-end');
addpath(toolboxPath,'-end');
pathdef_local = fullfile(userpath,'pathdef.m');
if isAdmin
    % Delete local pathdef file
    if isfile(pathdef_local)
        delete(pathdef_local);
    end
    % Save administrator local pathdef file
    savepath;
else
    % Create local user pathdef file
    fprintf('Updating local user "pathdef.m"...')
    savepath( pathdef_local );
    fprintf('[Complete]\n');
end

%% Rehash toolbox cache
fprintf('Rehashing Toolbox Cache...');
rehash TOOLBOXCACHE
fprintf('[Complete]\n');

%% internal functions (shared workspace)
% ------------------------------------------------------------------------
function migrateContent(sourceIn,destination,msg)

% Migrate toolbox folder contents
if ~isfolder(sourceIn)
    error(sprintf(...
        ['Change your working directory to the location of "install%s.m".\n',...
         '\n',...
         'If this problem persists:\n',...
         '\t(1) Unzip your original download of "%s" into a new directory\n',...
         '\t(2) Open a new instance of MATLAB "as administrator"\n',...
         '\t\t(a) Locate MATLAB shortcut\n',...
         '\t\t(b) Right click\n',...
         '\t\t(c) Select "Run as administrator"\n',...
         '\t(3) Change your "working directory" to the location of "install%s.m"\n',...
         '\t(4) Enter "install%s" (without quotes) into the command window\n',...
         '\t(5) Press Enter.'],toolboxShort,toolboxShort,toolboxShort,toolboxShort));
end

% Create Toolbox Path
[isDir,msgDir,msgID] = mkdir(destination);
if isDir
    fprintf('%s folder created successfully:\n\t"%s"\n',msg,destination);
else
    fprintf('Failed to create %s folder:\n\t"%s"\n',msg,destination);
    fprintf(adminSolution);
    error(msgID,msgDir);
end

% Migrate contents
files = dir(sourceIn);
wb = waitbar(0,sprintf('Copying %s contents...',msg));
n = numel(files);
fprintf('Copying %s contents:\n',msg);
for i = 1:n
    % source file location
    source = fullfile(sourceIn,files(i).name);
    
    if files(i).isdir
        switch files(i).name
            case '.'
                %Ignore
            case '..'
                %Ignore
            otherwise
                fprintf('\t%s...',files(i).name);
                nDestination = fullfile(destination,files(i).name);
                [isDir,msg,msgID] = mkdir(nDestination);
                if isDir
                    [isCopy,msg,msgID] = copyfile(source,nDestination,'f');
                    if isCopy
                        fprintf('[Complete]\n');
                    else
                        bin = msg == char(10);
                        msg(bin) = [];
                        bin = msg == char(13);
                        msg(bin) = [];
                        fprintf('[Failed: "%s"]\n',msg);
                    end
                else
                    bin = msg == char(10);
                    msg(bin) = [];
                    bin = msg == char(13);
                    msg(bin) = [];
                    fprintf('[Failed: "%s"]\n',msg);
                end
        end
    else
        fprintf('\t%s...',files(i).name);
        [isCopy,msg,msgID] = copyfile(source,destination,'f');
        
        if isCopy == 1
            fprintf('[Complete]\n');
        else
            bin = msg == char(10);
            msg(bin) = [];
            bin = msg == char(13);
            msg(bin) = [];
            fprintf('[Failed: "%s"]\n',msg);
        end
    end
    waitbar(i/n,wb);
end
set(wb,'Visible','off');
delete(wb);

end

end

%% Internal functions (unique workspace)
% ------------------------------------------------------------------------
function tfWrite = checkWriteAccess(pname)

tmpFname = fullfile(pname,'tmp.txt');
tmpHndle = fopen(tmpFname, 'w');
if tmpHndle < 0
    tfWrite = false;
else
    tfWrite = true;
    fclose(tmpHndle);
    delete(tmpFname);
end

end
% ------------------------------------------------------------------------
function removePath(toolboxName,pName,inPath,isPath,isDelete)
% Remove path
if inPath
    rmpath(pName);
    fprintf('%s path removed successfully:\n\t"%s"\n',toolboxName,pName);
end
% Remove folder
if isPath && isDelete
    [isRemoved, msg, msgID] = rmdir(pName,'s');
    if isRemoved
        fprintf('Previous version of %s removed successfully:\n\t"%s"\n',toolboxName,pName);
    else
        fprintf('Failed to remove old %s folder:\n\t"%s"\n',toolboxName,pName);
        %fprintf(adminSolution);
        error(msgID,msg);
    end
elseif ~isDelete
    fprintf('Skipping removal of old %s folder:\n\t"%s"\n',toolboxName,pName);
end

end