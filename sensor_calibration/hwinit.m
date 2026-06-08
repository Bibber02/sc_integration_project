scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = scriptFolder;
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(fullfile(projectRoot, 'hardware', 'rotating_pendulum'));

if exist('h', 'var') == 1 && isnumeric(h) && isscalar(h) && ~isempty(h)
    sampleTime = h;
else
    sampleTime = 0.01;
end

cfg = rotpend_hwinit(sampleTime);
h = cfg.h;
daoutoffs = cfg.daoutoffs;
daoutgain = cfg.daoutgain;
adinoffs = cfg.adinoffs;
adingain = cfg.adingain;
