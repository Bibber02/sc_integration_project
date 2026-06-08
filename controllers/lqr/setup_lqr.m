%% LQR setup
% Edit this block, then click Run.
%
% Angle convention used by this setup:
%   all-down: [pi; 0]
%   down-up:  [pi; pi]
%   all-up:   [0; 0]

targetEquilibrium = 'down_up';

switch targetEquilibrium
    case 'all_down'
        x0 = [pi; 0; 0; 0];
    case 'down_up'
        x0 = [pi; pi; 0; 0];
    case 'all_up'
        x0 = [0; 0; 0; 0];
    otherwise
        error('Unknown targetEquilibrium: %s', targetEquilibrium);
end

sampleTime = 0.001;

Q_lqr = diag([5 1 0.1 0.1]);
R_lqr = 5;

% The system-identification scripts use inputSign = -1, so the linear model
% input is the negative of the command signal sent to the real-time I/O block.
hardwareCommandToModelInputSign = -1;

openModelAfterSetup = true;

%% Setup logic
scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = scriptFolder;
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(projectRoot);
addpath(scriptFolder);

startupFile = fullfile(projectRoot, 'startup.m');
if isfile(startupFile)
    run(startupFile);
else
    scip.setupPath;
end
projectPaths = scip.paths;

settings = struct();
settings.targetEquilibrium = targetEquilibrium;
settings.x0 = x0;
settings.sampleTime = sampleTime;
settings.Q_lqr = Q_lqr;
settings.R_lqr = R_lqr;
settings.hardwareCommandToModelInputSign = hardwareCommandToModelInputSign;
settings.openModelAfterSetup = openModelAfterSetup;

lqrSettings = settings;

linearizedPlantFile = fullfile(projectPaths.lqr, 'lqr_linearized_plant.mat');
setupResultFile = fullfile(projectPaths.lqr, 'lqr_setup_result.mat');
simulinkModel = fullfile(projectPaths.lqr, 'rotating-pendulum', 'rotpentemplate.slx');
simulinkSupportFolder = fullfile(projectPaths.lqr, 'rotating-pendulum');

if isfolder(simulinkSupportFolder) && ~contains(path, simulinkSupportFolder)
    addpath(simulinkSupportFolder);
end

hiddenDefaults = struct();
hiddenDefaults.u0Mode = 'cancel_input_bias';
hiddenDefaults.hx = 1e-6;
hiddenDefaults.hu = 1e-6;
hiddenDefaults.filterOrder = 1;
hiddenDefaults.filterCutoffHz = 10;
hiddenDefaults.linearizedPlantFile = linearizedPlantFile;
hiddenDefaults.setupResultFile = setupResultFile;
hiddenDefaults.simulinkModel = simulinkModel;
hiddenDefaults.simulinkSupportFolder = simulinkSupportFolder;

linearizationSettings = settings;
linearizationSettings.u0Mode = hiddenDefaults.u0Mode;
linearizationSettings.hx = hiddenDefaults.hx;
linearizationSettings.hu = hiddenDefaults.hu;
linearizationSettings.saveOutput = true;
linearizationSettings.outputFile = linearizedPlantFile;

lin = linearize_rotpendulum(linearizationSettings);

A = lin.A;
B = lin.B;
C = lin.C;
D = lin.D;
Ad = lin.Ad;
Bd = lin.Bd;
Cd = lin.Cd;
Dd = lin.Dd;
sys_lin = lin.sys_lin;
sys_disc = lin.sys_disc;
Ts = lin.Ts;
x0 = lin.x0;
y0 = lin.y0;
u0 = lin.u0;
p = lin.p;

run(fullfile(scriptFolder, 'calc_lqr.m'));

if ~isscalar(hardwareCommandToModelInputSign) || hardwareCommandToModelInputSign == 0
    error('hardwareCommandToModelInputSign must be a nonzero scalar.');
end

K_lqr_command = -K_lqr / hardwareCommandToModelInputSign;
u_command0 = u0 / hardwareCommandToModelInputSign;

sys_kalman_disc = sys_disc;
sys_kalman_disc.B = hardwareCommandToModelInputSign * sys_disc.B;
sys_kalman_disc.D = hardwareCommandToModelInputSign * sys_disc.D;
sys_kalman_disc.InputName = {'u_command'};

h = sampleTime;
hardwareConfig = rotpend_hwinit(h);
daoutoffs = hardwareConfig.daoutoffs;
daoutgain = hardwareConfig.daoutgain;
adinoffs = hardwareConfig.adinoffs;
adingain = hardwareConfig.adingain;

[Q_kalman, R_kalman, P0_kalman, x0_kalman, kalmanObserverSettings] = ...
    loadKalmanObserverSettings(projectPaths, sampleTime);

setupResultFolder = fileparts(setupResultFile);
if ~isfolder(setupResultFolder)
    mkdir(setupResultFolder);
end

save(setupResultFile, ...
    'settings', 'hiddenDefaults', 'linearizationSettings', 'lin', ...
    'A', 'B', 'C', 'D', 'Ad', 'Bd', 'Cd', 'Dd', ...
    'sys_lin', 'sys_disc', 'sys_kalman_disc', 'Ts', 'x0', 'y0', 'u0', 'u_command0', 'p', ...
    'Q_lqr', 'R_lqr', 'K_lqr', 'K_lqr_command', 'closedLoopPoles', 'b', 'a', ...
    'hardwareCommandToModelInputSign', ...
    'Q_kalman', 'R_kalman', 'P0_kalman', 'x0_kalman', 'kalmanObserverSettings', ...
    'h', 'daoutoffs', 'daoutgain', 'adinoffs', 'adingain');

workspaceVars = {
    'settings'
    'A'
    'B'
    'C'
    'D'
    'Ad'
    'Bd'
    'Cd'
    'Dd'
    'sys_lin'
    'sys_disc'
    'sys_kalman_disc'
    'Ts'
    'x0'
    'y0'
    'u0'
    'u_command0'
    'p'
    'Q_lqr'
    'R_lqr'
    'K_lqr'
    'K_lqr_command'
    'closedLoopPoles'
    'hardwareCommandToModelInputSign'
    'b'
    'a'
    'Q_kalman'
    'R_kalman'
    'P0_kalman'
    'x0_kalman'
    'kalmanObserverSettings'
    'h'
    'daoutoffs'
    'daoutgain'
    'adinoffs'
    'adingain'
};

for k = 1:numel(workspaceVars)
    assignin('base', workspaceVars{k}, eval(workspaceVars{k}));
end

if openModelAfterSetup
    if ~isfile(simulinkModel)
        error('LQR Simulink model was not found: %s', simulinkModel);
    end

    if usejava('desktop')
        open_system(simulinkModel);
        fprintf('Opened LQR model: %s\n', simulinkModel);
    else
        load_system(simulinkModel);
        fprintf('Loaded LQR model: %s\n', simulinkModel);
    end

    configureKalmanFilterBlock(simulinkModel);
end

fprintf('\nLQR setup complete.\n');
fprintf('Linearized plant: %s\n', linearizedPlantFile);
fprintf('Setup result:     %s\n', setupResultFile);
fprintf('Target pose:      %s, x0 = [%g; %g; %g; %g]\n', ...
    targetEquilibrium, x0(1), x0(2), x0(3), x0(4));
fprintf('Sample time:      %.6g s\n', h);
fprintf('Command sign:     u_model = %.3g * u_command\n', hardwareCommandToModelInputSign);
fprintf('Kalman settings:  %s\n', kalmanObserverSettings.source);

function [Q_kalman, R_kalman, P0_kalman, x0_kalman, settings] = loadKalmanObserverSettings(projectPaths, sampleTime)
settings = struct();
settings.resultFile = fullfile(projectPaths.kalman, 'results', 'kalman_tuning_result.mat');
settings.noiseDataFile = fullfile(projectPaths.sensorNoiseMeasurement, 'noise_data.mat');
settings.source = 'default observer settings';
settings.sampleTime = sampleTime;

Q_kalman = diag([0.1, 0.1, 1, 0.1]);
R_kalman = defaultMeasurementCovariance(settings.noiseDataFile);
P0_kalman = 10 * eye(4);
x0_kalman = zeros(4, 1);

if ~isfile(settings.resultFile)
    return;
end

loaded = load(settings.resultFile);
requiredFields = {'Q_kalman', 'R_kalman', 'P0_kalman', 'x0_kalman'};
for kField = 1:numel(requiredFields)
    if ~isfield(loaded, requiredFields{kField})
        warning('Ignoring Kalman tuning result because it is missing %s: %s', ...
            requiredFields{kField}, settings.resultFile);
        return;
    end
end

Q_kalman = loaded.Q_kalman;
R_kalman = loaded.R_kalman;
P0_kalman = loaded.P0_kalman;
x0_kalman = loaded.x0_kalman;
settings.source = sprintf('Kalman tuning result (%s)', settings.resultFile);

if isfield(loaded, 'metrics') && isfield(loaded.metrics, 'selectedSource')
    settings.selectedSource = loaded.metrics.selectedSource;
end
end

function R = defaultMeasurementCovariance(noiseDataFile)
R = diag([7.3339511547e-6, 1.13475782135e-5]);

if ~isfile(noiseDataFile)
    return;
end

S = load(noiseDataFile, 'theta_1', 'theta_2');
if ~isfield(S, 'theta_1') || ~isfield(S, 'theta_2')
    return;
end

theta1 = localTimeseriesData(S.theta_1);
theta2 = localTimeseriesData(S.theta_2);
n = min(numel(theta1), numel(theta2));
noise = [theta1(1:n) - mean(theta1(1:n)), theta2(1:n) - mean(theta2(1:n))];
R = cov(noise, 1);
R = (R + R.') / 2;
end

function data = localTimeseriesData(value)
if isa(value, 'timeseries')
    data = value.Data;
else
    data = value;
end
data = double(squeeze(data));
data = data(:);
end

function configureKalmanFilterBlock(simulinkModel)
[~, modelName] = fileparts(simulinkModel);
kalmanBlock = [modelName, '/Kalman Filter'];
blockHandle = getSimulinkBlockHandle(kalmanBlock);

if blockHandle == -1
    warning('Kalman Filter block was not found in model %s.', modelName);
    return;
end

set_param(kalmanBlock, ...
    'ModelSourceVariable', 'sys_kalman_disc', ...
    'Q', 'Q_kalman', ...
    'R', 'R_kalman', ...
    'P0', 'P0_kalman', ...
    'X0', 'x0_kalman');

gainBlock = [modelName, '/Gain'];
if getSimulinkBlockHandle(gainBlock) ~= -1
    set_param(gainBlock, 'Gain', 'K_lqr_command');
else
    warning('LQR Gain block was not found in model %s.', modelName);
end

offsetBlock = [modelName, '/Constant'];
if getSimulinkBlockHandle(offsetBlock) ~= -1
    set_param(offsetBlock, 'Value', 'y0');
else
    warning('Output-offset Constant block was not found in model %s.', modelName);
end

fprintf('Configured in-memory Kalman Filter block to use Q_kalman/R_kalman and sys_kalman_disc.\n');
fprintf('Configured in-memory LQR Gain block to use K_lqr_command.\n');
fprintf('Configured in-memory output offset block to subtract y0.\n');
end
