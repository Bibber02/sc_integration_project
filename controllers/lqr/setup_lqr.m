%% Simple LQR + EKF setup
clear;
clc;

%% User settings
sampleTime = 0.01;

% Angle convention:
%   all-down: [pi; 0]
%   down-up:  [pi; pi]
%   all-up:   [0; 0]
x0 = [pi; pi; 0; 0];

Q_lqr = diag([5 1 0.1 0.1]);
R_lqr = 5;

%% Paths
scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptFolder));
modelFolder = fullfile(projectRoot, 'model');
hardwareFolder = fullfile(scriptFolder, 'rotating-pendulum');

addpath(modelFolder, '-begin');
addpath(scriptFolder, '-begin');
addpath(hardwareFolder, '-begin');

%% Hardware constants
run(fullfile(hardwareFolder, 'hwinit.m'));

%% Plant model and LQR
p = load_parameters();
lin = linearize_rotpendulum(struct( ...
    'x0', x0, ...
    'sampleTime', sampleTime, ...
    'p', p, ...
    'saveOutput', false));

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
h = Ts;
u0 = lin.u0;
y0 = lin.y0;
f0 = lin.f0;

[K_lqr, closedLoopPoles] = calc_lqr(sys_disc, Q_lqr, R_lqr);

%% EKF tuning result
[ekfResultFile, ekfTuning] = loadEkfTuning(projectRoot);
Q_ekf = ekfTuning.Q_ekf;
R_ekf = ekfTuning.R_ekf;
P0_ekf = ekfTuning.P0_ekf;
x0_ekf = ekfTuning.x0_ekf;
Ts_ekf = ekfTuning.Ts_ekf;

if abs(Ts - Ts_ekf) > 100 * eps(max(1, Ts_ekf))
    error('setup_lqr2:EkfSampleTimeMismatch', ...
        ['The LQR/EKF sample time is %.12g s, but the tuned EKF sample time is %.12g s. ', ...
        'Use sampleTime = Ts_ekf or retune the EKF.'], Ts, Ts_ekf);
end

% Compatibility aliases for old block/workspace names.
Q_kalman = Q_ekf;
R_kalman = R_ekf;
P0_kalman = P0_ekf;
x0_kalman = x0_ekf;

% Extra input passed to rotpendulumEkfStateTransition:
% [u_dev; Ts; u0; x0; p]
ekfInputParameters = [Ts; u0; x0; p];

% EKF block dialog values. The model folder was added to the path above, so
% Simulink can resolve these functions by name.
ekfStateTransitionFcn = 'rotpendulumEkfStateTransition';
ekfMeasurementFcn = 'rotpendulumEkfMeasurement';
ekfBlockSettings = struct( ...
    'StateTransitionFcn', ekfStateTransitionFcn, ...
    'StateTransitionHasExtraInput', true, ...
    'StateTransitionExtraInputSignal', '[u_dev; ekfInputParameters]', ...
    'MeasurementFcn', ekfMeasurementFcn, ...
    'ProcessNoise', 'Q_ekf', ...
    'MeasurementNoise', 'R_ekf', ...
    'InitialState', 'x0_ekf', ...
    'InitialStateCovariance', 'P0_ekf', ...
    'SampleTime', 'Ts');

fprintf('\nLQR/EKF setup complete.\n');
fprintf('Sample time:       %.6g s\n', Ts);
fprintf('EKF tuning result: %s\n', ekfResultFile);
fprintf('EKF state function: %s\n', ekfStateTransitionFcn);
fprintf('EKF measurement function: %s\n', ekfMeasurementFcn);

function [ekfResultFile, ekf] = loadEkfTuning(projectRoot)
candidates = {
    fullfile(projectRoot, 'kalman_filter_tuning', 'ekf_tuning_result.mat')
    fullfile(projectRoot, 'kalman_filter_tuning', 'EKF_tune_results.mat')
};

ekfResultFile = '';
for i = 1:numel(candidates)
    if isfile(candidates{i})
        ekfResultFile = candidates{i};
        break;
    end
end

if isempty(ekfResultFile)
    error('setup_lqr2:MissingEkfTuningResult', ...
        'Run kalman_filter_tuning/EKF_tune.m first.');
end

info = whos('-file', ekfResultFile);
fileNames = {info.name};
names = {'Q_ekf', 'R_ekf', 'P0_ekf', 'P0', 'x0_ekf', 'Ts_ekf', 'Ts'};
names = names(ismember(names, fileNames));
ekf = load(ekfResultFile, names{:});

if ~isfield(ekf, 'Q_ekf') || ~isequal(size(ekf.Q_ekf), [4 4])
    error('setup_lqr2:InvalidEkfQ', 'Q_ekf must be a 4-by-4 matrix.');
end
if ~isfield(ekf, 'R_ekf') || ~isequal(size(ekf.R_ekf), [2 2])
    error('setup_lqr2:InvalidEkfR', 'R_ekf must be a 2-by-2 matrix.');
end

if ~isfield(ekf, 'P0_ekf')
    if isfield(ekf, 'P0')
        ekf.P0_ekf = ekf.P0;
    else
        ekf.P0_ekf = 10 * eye(4);
    end
end
if ~isfield(ekf, 'x0_ekf')
    ekf.x0_ekf = zeros(4, 1);
end
if ~isfield(ekf, 'Ts_ekf')
    if isfield(ekf, 'Ts')
        ekf.Ts_ekf = ekf.Ts;
    else
        error('setup_lqr2:InvalidEkfTs', 'EKF result must contain Ts_ekf or Ts.');
    end
end

ekf.Q_ekf = double(ekf.Q_ekf);
ekf.R_ekf = double(ekf.R_ekf);
ekf.P0_ekf = double(ekf.P0_ekf);
ekf.x0_ekf = double(ekf.x0_ekf(:));
ekf.Ts_ekf = double(ekf.Ts_ekf);
end
