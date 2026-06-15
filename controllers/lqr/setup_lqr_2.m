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

%% LQR tuning settings
% State order:
%   x = [theta1; theta2; theta1_dot; theta2_dot]
%
% The LQR is designed for the discrete deviation model:
%   x_dev(k+1) = Ad*x_dev(k) + Bd*u_model_dev(k)
%
% Model-control law:
%   u_model_dev = -K_lqr*x_dev
%
% In Simulink, subtract the equilibrium from the EKF estimate first:
%   x_dev = x_hat - x_eq_lqr

useBrysonLqr = true;

% Manual fallback weights, used only when useBrysonLqr = false.
Q_lqr_manual = diag([5 3 0.1 0.01]);
R_lqr_manual = 1;

% Bryson-style allowed deviations. Smaller allowed deviation means larger
% LQR penalty for that state.
theta1_max = deg2rad(10);   % rad
theta2_max = deg2rad(8);    % rad
omega1_max = 2.0;           % rad/s
omega2_max = 2.0;           % rad/s

% Maximum desired control deviation. Increase this if the controller is too
% slow. Decrease this if the controller input is too aggressive.
maxControlDev = 0.40;

% Simple tuning multipliers. These are the main values to change.
q_theta1_mult = 1;
q_theta2_mult = 2;
q_omega1_mult = 1;
q_omega2_mult = 2;
r_mult        = 1;

% Sign convention between model input and actual command input.
%   u_model_dev = inputSignCommandToModel * u_command_dev
%
% Usually use +1. If your identified model/input convention is opposite to
% the real hardware command, set this to -1.
inputSignCommandToModel = 1;

% Optional linear closed-loop test before connecting Simulink.
runLqrDevSimulation  = true;
plotLqrDevSimulation = false;
N_lqrTest = 500;
xDev0_test = [deg2rad(5); deg2rad(5); 0; 0];
applyControlSaturationInTest = true;

if useBrysonLqr
    Q_lqr_base = diag(1 ./ [theta1_max theta2_max omega1_max omega2_max].^2);
    Q_lqr = diag([ ...
        q_theta1_mult * Q_lqr_base(1,1), ...
        q_theta2_mult * Q_lqr_base(2,2), ...
        q_omega1_mult * Q_lqr_base(3,3), ...
        q_omega2_mult * Q_lqr_base(4,4)]);
    R_lqr = r_mult * (1 / maxControlDev^2);
else
    Q_lqr = Q_lqr_manual;
    R_lqr = R_lqr_manual;
end

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

% Equilibrium values for Simulink.
x_eq_lqr = x0;
y_eq_lqr = y0;
u0_model = u0;
u0_command = inputSignCommandToModel * u0_model;

% LQR design for the model input convention.
[K_lqr, closedLoopPoles, lqrInfo] = calc_lqr(sys_disc, Q_lqr, R_lqr, Ts);

% Direct Simulink command gain.
%
% If your Simulink controller path computes command deviation directly from
% x_dev, use:
%   u_command_dev = K_lqr_command*x_dev
%
% If your Simulink plant/controller path uses model-input deviation, use:
%   u_model_dev = -K_lqr*x_dev
K_lqr_model = K_lqr;
K_lqr_command = -inputSignCommandToModel * K_lqr_model;

%% LQR diagnostics
n = size(Ad,1);
m = size(Bd,2);

openLoopPoles = eig(Ad);
closedLoopPoles = eig(Ad - Bd*K_lqr_model);
s_closedLoopPoles = log(closedLoopPoles)/Ts;
maxClosedLoopPoleMagnitude = max(abs(closedLoopPoles));
rankCtrb = rank(ctrb(Ad,Bd));

fprintf('\nLQR tuning settings:\n');
fprintf('  useBrysonLqr:               %d\n', useBrysonLqr);
fprintf('  theta1_max:                 %.6f rad = %.3f deg\n', theta1_max, rad2deg(theta1_max));
fprintf('  theta2_max:                 %.6f rad = %.3f deg\n', theta2_max, rad2deg(theta2_max));
fprintf('  omega1_max:                 %.6f rad/s\n', omega1_max);
fprintf('  omega2_max:                 %.6f rad/s\n', omega2_max);
fprintf('  maxControlDev:              %.6f\n', maxControlDev);
fprintf('  q multipliers:              [%.3g %.3g %.3g %.3g]\n', ...
    q_theta1_mult, q_theta2_mult, q_omega1_mult, q_omega2_mult);
fprintf('  r multiplier:               %.3g\n', r_mult);
fprintf('  inputSignCommandToModel:    %+d\n', inputSignCommandToModel);

fprintf('\nQ_lqr:\n');
disp(Q_lqr);
fprintf('R_lqr:\n');
disp(R_lqr);

fprintf('\nOpen-loop discrete poles:\n');
disp(openLoopPoles);

fprintf('LQR closed-loop poles in z-domain:\n');
disp(closedLoopPoles);

fprintf('LQR closed-loop poles in equivalent s-domain:\n');
disp(s_closedLoopPoles);

fprintf('Largest closed-loop pole magnitude: %.6f\n', maxClosedLoopPoleMagnitude);
fprintf('Controllability rank: %d / %d\n', rankCtrb, n);

if rankCtrb < n
    warning('The discrete plant is not fully controllable. LQR may still work if the system is stabilizable, but check the poles carefully.');
end

fprintf('\nK_lqr_model, for model law u_model_dev = -K_lqr_model*x_dev:\n');
disp(K_lqr_model);

fprintf('K_lqr_command, for direct Simulink command law u_command_dev = K_lqr_command*x_dev:\n');
disp(K_lqr_command);

%% Optional LQR-only simulation in deviation coordinates
if runLqrDevSimulation
    xLqrTest = zeros(n, N_lqrTest);
    uModelDevTest = zeros(m, N_lqrTest-1);
    uCommandDevTest = zeros(m, N_lqrTest-1);

    xLqrTest(:,1) = xDev0_test;

    for k = 1:N_lqrTest-1
        uModelDevTest(:,k) = -K_lqr_model*xLqrTest(:,k);
        uCommandDevTest(:,k) = inputSignCommandToModel * uModelDevTest(:,k);

        if applyControlSaturationInTest
            uCommandDevTest(:,k) = max(min(uCommandDevTest(:,k), maxControlDev), -maxControlDev);
            uModelDevTest(:,k) = inputSignCommandToModel * uCommandDevTest(:,k);
        end

        xLqrTest(:,k+1) = Ad*xLqrTest(:,k) + Bd*uModelDevTest(:,k);
    end

    fprintf('\nLQR-only deviation test:\n');
    fprintf('  Initial x_dev:              [% .4f % .4f % .4f % .4f]^T\n', xDev0_test);
    fprintf('  Final x_dev:                [% .4e % .4e % .4e % .4e]^T\n', xLqrTest(:,end));
    fprintf('  Max abs command deviation:  %.6f\n', max(abs(uCommandDevTest(:))));

    if plotLqrDevSimulation
        tLqrTest = (0:N_lqrTest-1)*Ts;

        figure;
        plot(tLqrTest, xLqrTest.');
        grid on;
        xlabel('Time [s]');
        ylabel('State deviation');
        legend('theta1','theta2','theta1 dot','theta2 dot');
        title('LQR-only linear deviation simulation');

        figure;
        stairs(tLqrTest(1:end-1), uCommandDevTest.');
        grid on;
        xlabel('Time [s]');
        ylabel('Command input deviation');
        title('LQR-only command deviation');
    end
end

%% Load EKF tuning result
ekfResultFile = fullfile(projectRoot, ...
    'kalman_filter_tuning', 'ekf_tuning_result.mat');

ekfTuning = load(ekfResultFile, ...
    'Q_ekf', 'R_ekf', 'P0_ekf', 'Ts_ekf');

Q_ekf = ekfTuning.Q_ekf;
R_ekf = ekfTuning.R_ekf;
P0_ekf = ekfTuning.P0_ekf;
Ts_ekf = ekfTuning.Ts_ekf;

if abs(Ts_ekf - Ts) > 1e-12
    warning('EKF tuning sample time Ts_ekf = %.16g differs from plant Ts = %.16g.', Ts_ekf, Ts);
end

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

% Additional parameters/input vector for the EKF block, if you use a single
% vector input for Ts and p.
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

fprintf('\nUse in Simulink LQR path:\n');
fprintf('  1. Compute x_dev = x_hat - x_eq_lqr.\n');
fprintf('  2. For model-input deviation, use u_model_dev = -K_lqr*x_dev.\n');
fprintf('  3. For direct command deviation, use u_command_dev = K_lqr_command*x_dev.\n');
fprintf('  4. Add u0_command if your plant/hardware expects absolute command input.\n');
