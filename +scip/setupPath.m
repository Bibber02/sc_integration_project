function setupPath()
%SETUPPATH Add shared project folders that do not create name conflicts.

p = scip.paths;
folders = {
    p.root
    p.lqr
    p.kalman
    p.hardwareRotatingPendulum
    p.model
    p.passiveLinkModels
    p.fullSystemGreyBoxViscous
    p.fullSystemGreyBoxStribeck
};

for k = 1:numel(folders)
    folder = folders{k};
    if isfolder(folder) && ~contains(path, folder)
        addpath(folder);
    end
end
end
