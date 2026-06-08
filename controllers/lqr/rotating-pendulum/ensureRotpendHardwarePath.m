function ensureRotpendHardwarePath(wrapperFile)
%ENSUREROTPENDHARDWAREPATH Add shared hardware and fugiboard support paths.

localFolder = fileparts(wrapperFile);
projectRoot = localFolder;
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end

addFolderIfMissing(localFolder, '-begin');
addFolderIfMissing(fullfile(projectRoot, 'hardware', 'rotating_pendulum'), '-end');

fugiboardFolders = {
    fullfile(localFolder, 'FUGIboardMatlab')
    fullfile(projectRoot, 'template', 'rotating-pendulum')
    fullfile(projectRoot, 'sensor_calibration')
    fullfile(projectRoot, 'system_identification', 'rotating-pendulum')
};

for k = 1:numel(fugiboardFolders)
    folder = fugiboardFolders{k};
    if isfolder(folder) && hasFugiboardMex(folder)
        addFolderIfMissing(folder, '-end');
        return;
    end
end
end

function tf = hasFugiboardMex(folder)
tf = ~isempty(dir(fullfile(folder, 'fugiboard.mex*'))) || ...
     isfile(fullfile(folder, 'fugiboard.m'));
end

function addFolderIfMissing(folder, position)
if ~isfolder(folder)
    return;
end

pathFolders = strsplit(path, pathsep);
if any(strcmpi(pathFolders, folder))
    return;
end

addpath(folder, position);
end
