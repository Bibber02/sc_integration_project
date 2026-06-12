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

Q_lqr = diag([5 3 0.1 0.01]);
R_lqr = 1;

%% Paths
scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptFolder));

modelFolder = fullfile(projectRoot, 'model');
hardwareFolder = fullfile(scriptFolder, 'rotating-pendulum');
ekfFolder = fullfile(projectRoot, 'extended_kalman_filter');

addpath(modelFolder, '-begin');
addpath(scriptFolder, '-begin');
addpath(hardwareFolder, '-begin');
addpath(ekfFolder, '-begin');

%% Check EKF functions
ekfStateFcnName = 'rotpendulumEkfStateTransition';
ekfMeasurementFcnName = 'rotpendulumEkfMeasurement';

if exist(ekfStateFcnName, 'file') ~= 2
    error('Could not find EKF state transition function: %s.m', ekfStateFcnName);
end

if exist(ekfMeasurementFcnName, 'file') ~= 2
    error('Could not find EKF measurement function: %s.m', ekfMeasurementFcnName);
end

ekfStateFcn = str2func(ekfStateFcnName);
ekfMeasurementFcn = str2func(ekfMeasurementFcnName);

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

%% Load EKF tuning result
ekfResultFile = fullfile(projectRoot, ...
    'kalman_filter_tuning', 'ekf_tuning_result.mat');

ekfTuning = load(ekfResultFile, ...
    'Q_ekf', 'R_ekf', 'P0_ekf', 'Ts_ekf');

Q_ekf = ekfTuning.Q_ekf;
R_ekf = ekfTuning.R_ekf;
P0_ekf = ekfTuning.P0_ekf;
Ts_ekf = ekfTuning.Ts_ekf;

%% EKF block variables

% Use these names in the Simulink EKF block:
%
%   State transition function: rotpendulumEkfStateTransition
%   Measurement function:      rotpendulumEkfMeasurement
%
% EKF settings:
%
%   Process noise:             Q_ekf
%   Measurement noise:         R_ekf
%   Initial state:             x0
%   Initial covariance:        P0_ekf
%   Sample time:               Ts
%
% Since rotpendulumEkfStateTransition has inputs:
%
%   x_k, u_k, Ts, p
%
% the EKF block supplies x_k internally.
% You must provide u_k, Ts, and p as additional inputs/parameters,
% depending on how your Simulink EKF block is configured.

ekfInputParameters = [Ts; p(:)];

fprintf('\nLQR/EKF setup complete.\n');
fprintf('Sample time:       %.6g s\n', Ts);
fprintf('EKF tuning result: %s\n', ekfResultFile);
fprintf('EKF folder:        %s\n', ekfFolder);
fprintf('State function:    %s\n', which(ekfStateFcnName));
fprintf('Measurement func:  %s\n', which(ekfMeasurementFcnName));
fprintf('\nUse in Simulink EKF block:\n');
fprintf('  State transition function: %s\n', ekfStateFcnName);
fprintf('  Measurement function:      %s\n', ekfMeasurementFcnName);