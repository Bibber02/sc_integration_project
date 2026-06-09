%% Single-file LQR setup
% Edit this block, then run this file from anywhere.

sampleTime = 0.01;

% Angle convention:
%   all-down: [pi; 0]
%   down-up:  [pi; pi]
%   all-up:   [0; 0]
x0 = [pi; pi; 0; 0];

Q_lqr = diag([5 1 0.1 0.1]);
R_lqr = 5;

filterOrder = 1;
filterCutoffHz = 10;

% Parameter loading / linearization settings.
l1 = 0.10;
g = 9.81;
hx = 1e-6;
hu = 1e-6;

% Hardware calibration constants used by the Simulink I/O subsystem.
daoutoffs = 0.00;
daoutgain = -6;
adinoffs = -[1.070535559513546 1.190099998448171];
adingain = [1.200374077783257 1.222725978474668];
adinoffs = [adinoffs 0 0 0 0 0];
adingain = [adingain 1 1 1 1 1];
sensorChannels = [6 7];
hardwareInitEnabled = true;

openModelAfterSetup = true;

%% Setup logic

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptFolder));
modelFolder = fullfile(projectRoot, 'model');
templateSupportFolder = fullfile(projectRoot, 'template', 'rotating-pendulum');
ekfResultCandidates = {
    fullfile(projectRoot, 'kalman_filter_tuning', 'ekf_tuning_result.mat')
    fullfile(projectRoot, 'kalman_filter_tuning', 'EKF_tune_results.mat')
};
ekfResultFile = selectEkfResultFile(ekfResultCandidates);
parameterResultFile = fullfile(projectRoot, ...
    'system_identification', 'full_system', 'grey_box', 'stribeck', ...
    'results_with_p0', 'full_system_id_new_data_locked_passive_result.mat');
linearizedPlantFile = fullfile(scriptFolder, 'lqr_linearized_plant.mat');
setupResultFile = fullfile(scriptFolder, 'lqr_setup_result.mat');
simulinkModel = fullfile(scriptFolder, 'rotating-pendulum', 'rotpentemplate.slx');

cd(projectRoot);
addpath(scriptFolder, '-begin');
addpath(modelFolder, '-begin');

ekfTuning = loadEkfTuningResult(ekfResultFile);
Q_ekf = ekfTuning.Q_ekf;
R_ekf = ekfTuning.R_ekf;
P0_ekf = ekfTuning.P0_ekf;
x0_ekf = ekfTuning.x0_ekf;
Ts_ekf = ekfTuning.Ts_ekf;

if abs(sampleTime - Ts_ekf) > 100 * eps(max(1, Ts_ekf))
    error('setup_lqr:EkfSampleTimeMismatch', ...
        ['The LQR/EKF sample time is %.12g s, but the tuned EKF sample time is %.12g s. ', ...
        'Run with sampleTime = Ts_ekf or retune the EKF at the desired sample time.'], ...
        sampleTime, Ts_ekf);
end

Q_kalman = Q_ekf;
R_kalman = R_ekf;
P0_kalman = P0_ekf;
x0_kalman = x0_ekf;

settings = struct();
settings.sampleTime = sampleTime;
settings.x0 = x0;
settings.Q_lqr = Q_lqr;
settings.R_lqr = R_lqr;
settings.filterOrder = filterOrder;
settings.filterCutoffHz = filterCutoffHz;
settings.l1 = l1;
settings.g = g;
settings.hx = hx;
settings.hu = hu;
settings.parameterResultFile = parameterResultFile;
settings.linearizedPlantFile = linearizedPlantFile;
settings.setupResultFile = setupResultFile;
settings.simulinkModel = simulinkModel;
settings.ekfResultFile = ekfResultFile;
settings.hardwareInitEnabled = hardwareInitEnabled;

lqrSettings = settings;
hiddenDefaults = struct();
hiddenDefaults.filterOrder = filterOrder;
hiddenDefaults.filterCutoffHz = filterCutoffHz;

p = loadPlantParameters(parameterResultFile, l1, g);
u0 = -p(6) / p(5);

[A, B, C, D, y0, f0] = linearizePlant(x0, u0, p, hx, hu);

sys_lin = ss(A, B, C, D);
sys_lin.StateName = {'theta1', 'theta2', 'theta1_dot', 'theta2_dot'};
sys_lin.InputName = {'u'};
sys_lin.OutputName = {'theta1', 'theta2'};

sys_disc = c2d(sys_lin, sampleTime, 'zoh');
Ad = sys_disc.A;
Bd = sys_disc.B;
Cd = sys_disc.C;
Dd = sys_disc.D;
Ts = sys_disc.Ts;
h = Ts;
ekfInputParameters = [Ts; u0; x0; p];
runHardwareInit = false; %#ok<NASGU>

run(fullfile(scriptFolder, 'calc_lqr.m'));

lin = struct();
lin.A = A;
lin.B = B;
lin.C = C;
lin.D = D;
lin.Ad = Ad;
lin.Bd = Bd;
lin.Cd = Cd;
lin.Dd = Dd;
lin.sys_lin = sys_lin;
lin.sys_disc = sys_disc;
lin.Ts = Ts;
lin.x0 = x0;
lin.u0 = u0;
lin.p = p;
lin.f0 = f0;
lin.y0 = y0;
lin.settings = settings;
lin.outputFile = linearizedPlantFile;

save(linearizedPlantFile, ...
    'A', 'B', 'C', 'D', 'Ad', 'Bd', 'Cd', 'Dd', ...
    'sys_lin', 'sys_disc', 'Ts', 'x0', 'u0', 'p', 'settings');

save(setupResultFile, ...
    'settings', 'hiddenDefaults', 'lin', ...
    'A', 'B', 'C', 'D', 'Ad', 'Bd', 'Cd', 'Dd', ...
    'sys_lin', 'sys_disc', 'Ts', 'x0', 'y0', 'u0', 'p', ...
    'Q_lqr', 'R_lqr', 'K_lqr', 'closedLoopPoles', 'b', 'a', ...
    'Q_ekf', 'R_ekf', 'P0_ekf', 'x0_ekf', 'Ts_ekf', ...
    'ekfInputParameters', 'ekfResultFile', ...
    'Q_kalman', 'R_kalman', 'P0_kalman', 'x0_kalman', ...
    'h', 'daoutoffs', 'daoutgain', 'adinoffs', 'adingain', ...
    'sensorChannels', 'hardwareInitEnabled');

workspaceVars = {
    'settings'
    'lqrSettings'
    'hiddenDefaults'
    'lin'
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
    'Ts'
    'h'
    'x0'
    'y0'
    'u0'
    'p'
    'Q_ekf'
    'R_ekf'
    'P0_ekf'
    'x0_ekf'
    'Ts_ekf'
    'ekfInputParameters'
    'ekfResultFile'
    'Q_lqr'
    'R_lqr'
    'K_lqr'
    'closedLoopPoles'
    'b'
    'a'
    'Q_kalman'
    'R_kalman'
    'P0_kalman'
    'x0_kalman'
    'daoutoffs'
    'daoutgain'
    'adinoffs'
    'adingain'
    'sensorChannels'
    'hardwareInitEnabled'
    'runHardwareInit'
};

for k = 1:numel(workspaceVars)
    assignin('base', workspaceVars{k}, eval(workspaceVars{k}));
end

load_system(simulinkModel);
[~, modelName] = fileparts(simulinkModel);
addpath(templateSupportFolder, '-begin');
configureSimulinkModel(modelName, sensorChannels);

if openModelAfterSetup && usejava('desktop')
    open_system(modelName);
else
    close_system(modelName, 0);
end

runHardwareInit = hardwareInitEnabled;
assignin('base', 'runHardwareInit', runHardwareInit);

fprintf('\nLQR setup complete.\n');
fprintf('Project root:      %s\n', projectRoot);
fprintf('Sample time:       %.6g s\n', Ts);
fprintf('Initial state x0:  [%g; %g; %g; %g]\n', x0(1), x0(2), x0(3), x0(4));
fprintf('Linearized plant:  %s\n', linearizedPlantFile);
fprintf('Setup result:      %s\n', setupResultFile);
fprintf('Simulink model:    %s\n', simulinkModel);

function ekfResultFile = selectEkfResultFile(ekfResultCandidates)
ekfResultFile = ekfResultCandidates{1};
for i = 1:numel(ekfResultCandidates)
    if isfile(ekfResultCandidates{i})
        ekfResultFile = ekfResultCandidates{i};
        return;
    end
end
end

function ekfTuning = loadEkfTuningResult(ekfResultFile)
if ~isfile(ekfResultFile)
    error('setup_lqr:MissingEkfTuningResult', ...
        ['Could not find the EKF tuning result file:\n  %s\n', ...
        'Run kalman_filter_tuning/EKF_tune.m first.'], ekfResultFile);
end

namesToLoad = variablesToLoad(ekfResultFile);
ekfTuning = load(ekfResultFile, namesToLoad{:});
ekfTuning = normalizeEkfTuningResult(ekfTuning, ekfResultFile);

ekfTuning.Q_ekf = double(ekfTuning.Q_ekf);
ekfTuning.R_ekf = double(ekfTuning.R_ekf);
ekfTuning.P0_ekf = double(ekfTuning.P0_ekf);
ekfTuning.x0_ekf = double(ekfTuning.x0_ekf(:));
ekfTuning.Ts_ekf = double(ekfTuning.Ts_ekf);

if ~isequal(size(ekfTuning.Q_ekf), [4 4])
    error('setup_lqr:InvalidEkfQ', ...
        'Q_ekf must be a 4-by-4 matrix in %s.', ekfResultFile);
end
if ~isequal(size(ekfTuning.R_ekf), [2 2])
    error('setup_lqr:InvalidEkfR', ...
        'R_ekf must be a 2-by-2 matrix in %s.', ekfResultFile);
end
if ~isequal(size(ekfTuning.P0_ekf), [4 4])
    error('setup_lqr:InvalidEkfP0', ...
        'P0_ekf must be a 4-by-4 matrix in %s.', ekfResultFile);
end
if numel(ekfTuning.x0_ekf) ~= 4
    error('setup_lqr:InvalidEkfX0', ...
        'x0_ekf must contain four states in %s.', ekfResultFile);
end
if ~isscalar(ekfTuning.Ts_ekf) || ~isfinite(ekfTuning.Ts_ekf) || ekfTuning.Ts_ekf <= 0
    error('setup_lqr:InvalidEkfTs', ...
        'Ts_ekf must be a positive scalar in %s.', ekfResultFile);
end
end

function names = variablesToLoad(ekfResultFile)
fileInfo = whos('-file', ekfResultFile);
fileNames = {fileInfo.name};
candidateNames = {'Q_ekf', 'R_ekf', 'P0_ekf', 'P0', 'x0_ekf', 'Ts_ekf', 'Ts'};
names = candidateNames(ismember(candidateNames, fileNames));
end

function ekfTuning = normalizeEkfTuningResult(ekfTuning, ekfResultFile)
if ~isfield(ekfTuning, 'Q_ekf')
    error('setup_lqr:InvalidEkfTuningResult', ...
        'EKF tuning result %s is missing variable Q_ekf.', ekfResultFile);
end
if ~isfield(ekfTuning, 'R_ekf')
    error('setup_lqr:InvalidEkfTuningResult', ...
        'EKF tuning result %s is missing variable R_ekf.', ekfResultFile);
end
if ~isfield(ekfTuning, 'P0_ekf')
    if isfield(ekfTuning, 'P0')
        ekfTuning.P0_ekf = ekfTuning.P0;
    else
        ekfTuning.P0_ekf = 10 * eye(4);
    end
end
if ~isfield(ekfTuning, 'x0_ekf')
    ekfTuning.x0_ekf = zeros(4, 1);
end
if ~isfield(ekfTuning, 'Ts_ekf')
    if isfield(ekfTuning, 'Ts')
        ekfTuning.Ts_ekf = ekfTuning.Ts;
    else
        error('setup_lqr:InvalidEkfTuningResult', ...
            'EKF tuning result %s is missing variable Ts_ekf or Ts.', ...
            ekfResultFile);
    end
end
end

function p = loadPlantParameters(parameterResultFile, l1, g)
S = load(parameterResultFile, 'resultTable');
values = S.resultTable{:, 2};

if iscell(values)
    values = cellfun(@double, values);
elseif isstring(values) || ischar(values)
    values = str2double(values);
end

values = double(values(:));

pa       = values(1);
pb1      = values(2);
pc1      = values(3);
pg1      = values(4);
pu       = values(5);
p0       = 0;
pb2      = values(7);
pg2      = values(8);
pc2      = values(9);
psdelta2 = values(10);
vs2      = values(11);
epsv1    = values(12);
epsv2    = values(13);
pc       = (l1 / g) * pg2;

p = [
    pa
    pb1
    pc1
    pg1
    pu
    p0
    pb2
    pg2
    pc2
    psdelta2
    vs2
    epsv1
    epsv2
    pc
];
end

function [A, B, C, D, y0, f0] = linearizePlant(x0, u0, p, hx, hu)
[f0, y0] = nonlinearPlant(x0, u0, p);

nx = numel(x0);
nu = 1;
ny = numel(y0);

A = zeros(nx, nx);
B = zeros(nx, nu);
C = zeros(ny, nx);
D = zeros(ny, nu);

for i = 1:nx
    dx = zeros(nx, 1);
    dx(i) = hx;

    [fPlus, yPlus] = nonlinearPlant(x0 + dx, u0, p);
    [fMinus, yMinus] = nonlinearPlant(x0 - dx, u0, p);

    A(:, i) = (fPlus - fMinus) / (2 * hx);
    C(:, i) = (yPlus - yMinus) / (2 * hx);
end

[fPlus, yPlus] = nonlinearPlant(x0, u0 + hu, p);
[fMinus, yMinus] = nonlinearPlant(x0, u0 - hu, p);

B(:, 1) = (fPlus - fMinus) / (2 * hu);
D(:, 1) = (yPlus - yMinus) / (2 * hu);
end

function [xdot, y] = nonlinearPlant(x, u, p)
theta1 = x(1);
theta2 = x(2);
dtheta1 = x(3);
dtheta2 = x(4);

pa = p(1);
pb1 = p(2);
pc1 = p(3);
pg1 = p(4);
pu = p(5);
p0 = p(6);
pb2 = p(7);
pg2 = p(8);
pc2 = p(9);
psdelta2 = p(10);
vs2 = p(11);
epsv1 = p(12);
epsv2 = p(13);
pc = p(14);

c2 = cos(theta2);
s2 = sin(theta2);

M11 = pa + 1 + 2 * pc * c2;
M12 = 1 + pc * c2;
M21 = M12;
M22 = 1;

Cqdot1 = -pc * s2 * (2 * dtheta1 * dtheta2 + dtheta2^2);
Cqdot2 = pc * s2 * dtheta1^2;

F1 = pb1 * dtheta1 + pc1 * tanh(dtheta1 / epsv1);

stribeck2 = pc2 + psdelta2 * exp(-(dtheta2 / vs2)^2);
F2 = pb2 * dtheta2 + stribeck2 * tanh(dtheta2 / epsv2);

G1 = -pg1 * sin(theta1) - pg2 * sin(theta1 + theta2);
G2 = -pg2 * sin(theta1 + theta2);

I1 = pu * u + p0;
I2 = 0;

rhs1 = I1 - Cqdot1 - F1 - G1;
rhs2 = I2 - Cqdot2 - F2 - G2;

detM = M11 * M22 - M12 * M21;

ddtheta1 = (M22 * rhs1 - M12 * rhs2) / detM;
ddtheta2 = (-M21 * rhs1 + M11 * rhs2) / detM;

xdot = zeros(4, 1);
xdot(1) = dtheta1;
xdot(2) = dtheta2;
xdot(3) = ddtheta1;
xdot(4) = ddtheta2;

y = zeros(2, 1);
y(1) = theta1;
y(2) = theta2;
end

function configureSimulinkModel(modelName, sensorChannels)
set_param([modelName '/Gain'], 'Gain', '-K_lqr');
set_param([modelName '/Constant'], 'Value', 'y0');

configureExtendedKalmanFilter(modelName);

set_param([modelName '/Subsystem/RT input'], ...
    'Ts', 'h', ...
    'channels', mat2str(sensorChannels), ...
    'MaskInitialization', hardwareMaskInitialization());
end

function configureExtendedKalmanFilter(modelName)
deleteBlockIfPresent([modelName '/Kalman Filter']);
deleteBlockIfPresent([modelName '/Extended Kalman Filter']);
deleteBlockIfPresent([modelName '/EKF Input Vector']);
deleteBlockIfPresent([modelName '/EKF Parameters']);

load_system('ctrlSharedLib');

ekfBlock = [modelName '/Extended Kalman Filter'];
add_block('ctrlSharedLib/Extended Kalman Filter', ekfBlock, ...
    'Position', [535 204 625 266]);

set_param(ekfBlock, ...
    'StateTransitionFcn', 'rotpendulumEkfStateTransition', ...
    'HasStateTransitionJacobianFcn', 'off', ...
    'HasAdditiveProcessNoise', 'Additive', ...
    'ProcessNoise', 'Q_ekf', ...
    'HasTimeVaryingProcessNoise', 'off', ...
    'HasStateTransitionFcnExtraArgument', '1', ...
    'InitialState', 'x0_ekf', ...
    'InitialStateCovariance', 'P0_ekf', ...
    'MeasurementFcn1', 'rotpendulumEkfMeasurement', ...
    'HasMeasurementJacobianFcn1', 'off', ...
    'HasAdditiveMeasurementNoise1', 'Additive', ...
    'MeasurementNoise1', 'R_ekf', ...
    'HasTimeVaryingMeasurementNoise1', 'off', ...
    'HasMeasurementFcnExtraArgument1', '0', ...
    'EnableMultirate', 'off', ...
    'UseCurrentEstimator', 'on', ...
    'OutputStateCovariance', 'off', ...
    'SampleTime', 'Ts');

paramBlock = [modelName '/EKF Parameters'];
add_block('simulink/Sources/Constant', paramBlock, ...
    'Position', [330 340 455 370], ...
    'Value', 'ekfInputParameters');

muxBlock = [modelName '/EKF Input Vector'];
add_block('simulink/Signal Routing/Mux', muxBlock, ...
    'Position', [485 315 490 365], ...
    'Inputs', '2');

connectLine(modelName, 'Saturation/1', 'EKF Input Vector/1');
connectLine(modelName, 'EKF Parameters/1', 'EKF Input Vector/2');
connectLine(modelName, 'Sum1/1', 'Extended Kalman Filter/1');
connectLine(modelName, 'EKF Input Vector/1', 'Extended Kalman Filter/2');
connectLine(modelName, 'Extended Kalman Filter/1', 'Gain/1');
connectLine(modelName, 'Extended Kalman Filter/1', 'Demux1/1');
end

function deleteBlockIfPresent(blockPath)
if blockExists(blockPath)
    delete_block(blockPath);
end
end

function tf = blockExists(blockPath)
tf = getSimulinkBlockHandle(blockPath) ~= -1;
end

function connectLine(modelName, sourcePort, destinationPort)
deleteLineAtDestination(modelName, destinationPort);
add_line(modelName, sourcePort, destinationPort, 'autorouting', 'on');
end

function deleteLineAtDestination(modelName, destinationPort)
[blockName, portNumber] = splitBlockPort(destinationPort);
blockPath = [modelName '/' blockName];
portHandles = get_param(blockPath, 'PortHandles');
lineHandle = get_param(portHandles.Inport(portNumber), 'Line');
if lineHandle ~= -1
    delete_line(lineHandle);
end
end

function [blockName, portNumber] = splitBlockPort(blockPort)
slashIndex = find(blockPort == '/', 1, 'last');
blockName = blockPort(1:slashIndex - 1);
portNumber = str2double(blockPort(slashIndex + 1:end));
end

function code = hardwareMaskInitialization()
code = sprintf([ ...
    'Ts = evalin(''base'', ''Ts'');\n' ...
    'h = evalin(''base'', ''h'');\n' ...
    'daoutoffs = evalin(''base'', ''daoutoffs'');\n' ...
    'daoutgain = evalin(''base'', ''daoutgain'');\n' ...
    'adinoffs = evalin(''base'', ''adinoffs'');\n' ...
    'adingain = evalin(''base'', ''adingain'');\n' ...
    'runHardwareInit = evalin(''base'', ''exist(''''runHardwareInit'''', ''''var'''') && runHardwareInit'');\n' ...
    'if runHardwareInit\n' ...
    'fugihandle = fugiboard(''Open'', ''Pendulum1'');\n' ...
    'fugihandle.WatchdogTimeout = 0.5;\n' ...
    'fugiboard(''SetParams'', fugihandle);\n' ...
    'fugiboard(''Write'', fugihandle, 0, 0, 0, 0);\n' ...
    'fugiboard(''Write'', fugihandle, 5, 1, 0, 0);\n' ...
    'data = fugiboard(''Read'', fugihandle);\n' ...
    'fpgaModel = bitshift(data(1), -5);\n' ...
    'fpgaVersion = bitand(data(1), 31);\n' ...
    'disp(sprintf(''FPGA setup %%d, version %%d'', fpgaModel, fpgaVersion));\n' ...
    'fugiboard(''Write'', fugihandle, 0, 1, 0, 0);\n' ...
    'pause(0.1);\n' ...
    'else\n' ...
    'fugihandle = [];\n' ...
    'end\n']);
end
