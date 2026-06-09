function [sys,x0,str,ts] = sfusbout(t,x,u,flag,Ts,h)
ensureRotpendSharedPath(mfilename('fullpath'));
[sys,x0,str,ts] = rotpend_sfusbout(t,x,u,flag,Ts,h);
end

function ensureRotpendSharedPath(wrapperFile)
projectRoot = fileparts(wrapperFile);
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(fullfile(projectRoot, 'hardware', 'rotating_pendulum'));
end
