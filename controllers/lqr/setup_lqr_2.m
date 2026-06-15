%% Simple LQR + EKF setup
clear;
clc;

%% User settings
Ts = 0.01;
sampleTime = Ts;

% Equilibrium / linearization point.
% Angle convention:
%   all-down: [pi; 0]
%   down-up:  [pi; pi]
%   all-up:   [0; 0]
x0 = [pi; pi; 0; 0];

%% LQR tuning
% State order:
%   x = [theta1; theta2; theta1_dot; theta2_dot]
%
% MAIN TUNING RULES:
%   Increase R_manual  -> less aggressive input, smaller K_lqr
%   Decrease R_manual  -> more aggressive input, larger K_lqr
%   Increase Q(2,2)    -> care more about theta2 angle
%   Increase Q(4,4)    -> care more about theta2 velocity / damping
%
% Start simple. Tune mostly R_manual first.
Q_manual = diag([5 3 0.1 0.01]);
R_manual = 100;

% Use this only if the motor command sign is opposite to the model input sign.
% Usually keep this at +1.
inputSignCommandToModel = 1;

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

%% EKF function names for Simulink
ekfStateFcnName = 'rotpendulumEkfStateTransition';
ekfMeasurementFcnName = 'rotpendulumEkfMeasurement';

if exist(ekfStateFcnName, 'file') ~= 2
    error('Could not find %s.m', ekfStateFcnName);
end

if exist(ekfMeasurementFcnName, 'file') ~= 2
    error('Could not find %s.m', ekfMeasurementFcnName);
end

%% Hardware constants
run(fullfile(hardwareFolder, 'hwinit.m'));

%% Linearized plant model
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

%% LQR design
Q_lqr = Q_manual;
R_lqr = R_manual;

[K_lqr, ~, closedLoopPoles] = dlqr(Ad, Bd, Q_lqr, R_lqr);

Acl = Ad - Bd*K_lqr;
openLoopPoles = eig(Ad);
closedLoopPoles = eig(Acl);

openLoopPoles_s = log(openLoopPoles)/Ts;
closedLoopPoles_s = log(closedLoopPoles)/Ts;

% Gain to use directly in a Simulink Gain block if the block input is x_dev.
% This gives:
%   u_command_dev = K_lqr_simulink*x_dev
%
% With inputSignCommandToModel = +1, this is simply -K_lqr.
K_lqr_simulink = -inputSignCommandToModel*K_lqr;

% Equilibrium values for Simulink.
% Use:
%   x_dev = x_hat - x_eq_lqr
%   u_command = u0_command + K_lqr_simulink*x_dev
x_eq_lqr = x0;
y_eq_lqr = y0;
u0_model = u0;
u0_command = inputSignCommandToModel*u0_model;

%% Load EKF tuning result
ekfResultFile = fullfile(projectRoot, ...
    'extended_kalman_filter', 'ekf_tuning_result.mat');

ekfTuning = load(ekfResultFile, ...
    'Q_ekf', 'R_ekf', 'P0_ekf', 'Ts_ekf');

Q_ekf = ekfTuning.Q_ekf;
R_ekf = ekfTuning.R_ekf;
P0_ekf = ekfTuning.P0_ekf;
Ts_ekf = ekfTuning.Ts_ekf;

% Parameters for the EKF block.
ekfInputParameters = [Ts; p(:)];

%% Minimal diagnostics
[~, dominantIdx] = max(abs(closedLoopPoles));

disp('================ LQR SETUP ================');
fprintf('Ts = %.6f s\n', Ts);
fprintf('R_lqr = %.6g\n', R_lqr);
fprintf('Q_lqr diagonal = ');
disp(diag(Q_lqr).');

fprintf('Controllability rank = %d / %d\n', rank(ctrb(Ad,Bd)), size(Ad,1));

fprintf('\nK_lqr for model law u_model_dev = -K_lqr*x_dev:\n');
disp(K_lqr);

fprintf('K_lqr_simulink for direct Gain block u_command_dev = K_lqr_simulink*x_dev:\n');
disp(K_lqr_simulink);

fprintf('\nOpen-loop poles in z-domain:\n');
disp(openLoopPoles);

fprintf('Closed-loop poles in z-domain:\n');
disp(closedLoopPoles);

fprintf('Closed-loop equivalent poles in s-domain:\n');
disp(closedLoopPoles_s);

fprintf('Dominant closed-loop pole: z = %.6g%+.6gi, s = %.6g%+.6gi\n', ...
    real(closedLoopPoles(dominantIdx)), imag(closedLoopPoles(dominantIdx)), ...
    real(closedLoopPoles_s(dominantIdx)), imag(closedLoopPoles_s(dominantIdx)));

fprintf('\nIf the controller is too aggressive, increase R_manual.\n');
fprintf('If theta2 is too oscillatory, increase Q_manual(2,2) or Q_manual(4,4).\n');
disp('===========================================');
