function observerModel = kalman_load_linear_observer_model(projectPaths, settings)
%KALMAN_LOAD_LINEAR_OBSERVER_MODEL Load the linear model used by the observer.

candidates = {
    fullfile(projectPaths.lqr, 'lqr_setup_result.mat')
    fullfile(projectPaths.lqr, 'lqr_linearized_plant.mat')
    projectPaths.linearizedPlantUprightTs001
};

modelData = [];
sourceFile = '';
for k = 1:numel(candidates)
    if isfile(candidates{k})
        modelData = load(candidates{k});
        sourceFile = candidates{k};
        break;
    end
end

if isempty(modelData)
    linSettings = struct();
    % Angle convention: all-down [pi;0], down-up [pi;pi], all-up [0;0].
    linSettings.x0 = [pi; pi; 0; 0];
    linSettings.sampleTime = settings.runtimeSampleTime;
    linSettings.saveOutput = false;
    lin = linearize_rotpendulum(linSettings);
    modelData = lin;
    sourceFile = 'linearize_rotpendulum';
end

if isfield(modelData, 'sys_lin')
    sysData = c2d(modelData.sys_lin, settings.TsData, 'zoh');
    runtimeSampleTime = getRuntimeSampleTime(modelData, settings);
    sysRuntime = c2d(modelData.sys_lin, runtimeSampleTime, 'zoh');
elseif isfield(modelData, 'sys_disc')
    sysData = d2d(modelData.sys_disc, settings.TsData, 'zoh');
    runtimeSampleTime = getRuntimeSampleTime(modelData, settings);
    sysRuntime = d2d(modelData.sys_disc, runtimeSampleTime, 'zoh');
else
    error('Could not find sys_lin or sys_disc in the observer model source.');
end

if isfield(modelData, 'x0')
    x0 = modelData.x0(:);
elseif isfield(modelData, 'lin') && isfield(modelData.lin, 'x0')
    x0 = modelData.lin.x0(:);
else
    x0 = [pi; pi; 0; 0];
end

if isfield(modelData, 'u0')
    u0 = modelData.u0;
elseif isfield(modelData, 'lin') && isfield(modelData.lin, 'u0')
    u0 = modelData.lin.u0;
else
    u0 = 0;
end
u0 = double(u0(1));

if isfield(modelData, 'lin') && isfield(modelData.lin, 'y0')
    y0 = modelData.lin.y0(:);
else
    y0 = sysData.C * x0 + sysData.D * u0;
end

observerModel = struct();
observerModel.sourceFile = sourceFile;
observerModel.sysData = sysData;
observerModel.sysRuntime = sysRuntime;
observerModel.runtimeSampleTime = runtimeSampleTime;
observerModel.x0 = x0;
observerModel.u0 = u0;
observerModel.y0 = y0(:);
end

function runtimeSampleTime = getRuntimeSampleTime(modelData, settings)
if isfield(modelData, 'Ts') && isfinite(modelData.Ts) && modelData.Ts > 0
    runtimeSampleTime = modelData.Ts;
elseif isfield(modelData, 'sys_disc') && isfinite(modelData.sys_disc.Ts) && modelData.sys_disc.Ts > 0
    runtimeSampleTime = modelData.sys_disc.Ts;
else
    runtimeSampleTime = settings.runtimeSampleTime;
end
end
