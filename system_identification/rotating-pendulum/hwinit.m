scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = scriptFolder;
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(fullfile(projectRoot, 'hardware', 'rotating_pendulum'));

cfg = rotpend_hwinit(0.01);
h = cfg.h;
daoutoffs = cfg.daoutoffs;
daoutgain = cfg.daoutgain;
adinoffs = cfg.adinoffs;
adingain = cfg.adingain;
