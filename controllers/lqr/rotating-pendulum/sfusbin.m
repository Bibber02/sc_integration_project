function [sys,x0,str,ts] = sfusbin(t,x,u,flag,Ts,h,channels)
ensureRotpendSharedPath(mfilename('fullpath'));
[sys,x0,str,ts] = rotpend_sfusbin(t,x,u,flag,Ts,h,channels);
end

function ensureRotpendSharedPath(wrapperFile)
projectRoot = fileparts(wrapperFile);
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(fullfile(projectRoot, 'hardware', 'rotating_pendulum'));
end
