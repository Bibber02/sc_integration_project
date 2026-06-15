%% Simple LQI + EKF setup
clear;
clc;

%% User settings
sampleTime = 0.01;

% Angle convention:
%   all-down: [pi; 0]
%   down-up:  [pi; pi]
%   all-up:   [0; 0]
x0 = [pi; pi; 0; 0];

%% LQI tuning settings
% State order:
%   x = [theta1; theta2; theta1_dot; theta2_dot]
%
% Deviation model:
%   x_dev(k+1) = Ad*x_dev(k) + Bd*u_model_dev(k)
%
% Integrator:
%   xi(k+1) = xi(k) + Ts*(r_dev(k) - y_track_dev(k))
%
% Model-input controller:
%   u_model_dev = -Kx_lqi*x_dev - Ki_lqi*xi + G_ref*r_dev

% Output to track with the integrator and reference feedforward.
% theta1 only:
C_track = [1 0 0 0];
D_track = 0;

% Main tuning knobs.
% Increase R_lqi if the controller is too aggressive.
% Increase Qi_lqi if the reference error disappears too slowly.
Qx_lqi = diag([100 100 100 0.1]);
Qi_lqi = 1;
R_lqi  = 1;

% Sign convention between model input and actual command input.
%   u_model_dev = inputSignCommandToModel * u_command_dev
% Usually use +1. Use -1 only if your command sign is opposite to the
% identified model input sign.
inputSignCommandToModel = 1;

% Optional safety value used only by the simple linear test below.
maxControlDev = 0.40;

% Optional linear LQI test in deviation coordinates.
runLqiLinearTest = false;
plotLqiLinearTest = false;
N_lqiTest = 500;
rDev_test = deg2rad(5);
xDev0_test = zeros(4,1);
xi0_test = 0;

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

%% EKF function names
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

%% Plant model
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

%% Equilibrium values for Simulink
x_eq_lqi = x0;
y_eq_lqi = y0;
u0_model = u0;
u0_command = inputSignCommandToModel * u0_model;

% Reference offset for the tracked output.
% For theta1 tracking, use:
%   r_dev = theta1_ref_abs - theta1_0
%   y_track_dev = theta1_measured_or_estimated - theta1_0
theta1_0 = C_track*x0 + D_track*u0;
y_track_0 = theta1_0;

%% LQI design
n = size(Ad,1);
m = size(Bd,2);
p_int = size(C_track,1);

A_aug = [Ad,              zeros(n,p_int);
        -Ts*C_track,      eye(p_int)];

B_aug = [Bd;
         zeros(p_int,m)];

Q_aug = blkdiag(Qx_lqi, Qi_lqi);
R_aug = R_lqi;

K_aug = dlqr(A_aug, B_aug, Q_aug, R_aug);

Kx_lqi = K_aug(:,1:n);
Ki_lqi = K_aug(:,n+1:end);

closedLoopPoles_lqi = eig(A_aug - B_aug*K_aug);
s_closedLoopPoles_lqi = log(closedLoopPoles_lqi)/Ts;

%% Reference feedforward matrix
% Regulator equations for steady-state reference tracking:
%   x_ss = Pi_ref*r_dev
%   u_ss = Gamma_ref*r_dev
%
% With the controller:
%   u_model_dev = -Kx_lqi*x_dev - Ki_lqi*xi + G_ref*r_dev
%
% this feedforward makes the reference response easier for the integrator.
% The integrator still removes remaining steady-state error.
M_ref = [Ad - eye(n), Bd;
         C_track,     D_track];

rhs_ref = [zeros(n,p_int);
           eye(p_int)];

X_ref = M_ref \ rhs_ref;

Pi_ref = X_ref(1:n,:);
Gamma_ref = X_ref(n+1:n+m,:);

G_ref_model = Gamma_ref + Kx_lqi*Pi_ref;
G_ref_command = inputSignCommandToModel * G_ref_model;

res_ref = M_ref*X_ref - rhs_ref;
rel_res_ref = norm(res_ref) / max(1,norm(rhs_ref));

%% Simulink controller gains
% Model-input form:
%   u_model_dev = -Kx_lqi*x_dev - Ki_lqi*xi + G_ref_model*r_dev
%
% Direct command-input form:
%   u_command_dev = Kx_lqi_command*x_dev + Ki_lqi_command*xi + G_ref_command*r_dev
Kx_lqi_model = Kx_lqi;
Ki_lqi_model = Ki_lqi;

Kx_lqi_command = -inputSignCommandToModel * Kx_lqi_model;
Ki_lqi_command = -inputSignCommandToModel * Ki_lqi_model;

K_lqi_model = [Kx_lqi_model Ki_lqi_model];
K_lqi_command = [Kx_lqi_command Ki_lqi_command];

%% Optional linear LQI test
if runLqiLinearTest
    xTest = zeros(n,N_lqiTest);
    xiTest = zeros(p_int,N_lqiTest);
    uModelDevTest = zeros(m,N_lqiTest-1);
    uCommandDevTest = zeros(m,N_lqiTest-1);

    xTest(:,1) = xDev0_test;
    xiTest(:,1) = xi0_test;

    for k = 1:N_lqiTest-1
        yTrackDev = C_track*xTest(:,k);
        eTrack = rDev_test - yTrackDev;

        uModelDevTest(:,k) = -Kx_lqi_model*xTest(:,k) ...
                             -Ki_lqi_model*xiTest(:,k) ...
                             +G_ref_model*rDev_test;

        uCommandDevTest(:,k) = inputSignCommandToModel*uModelDevTest(:,k);
        uCommandDevTest(:,k) = max(min(uCommandDevTest(:,k), maxControlDev), -maxControlDev);
        uModelDevTest(:,k) = inputSignCommandToModel*uCommandDevTest(:,k);

        xTest(:,k+1) = Ad*xTest(:,k) + Bd*uModelDevTest(:,k);
        xiTest(:,k+1) = xiTest(:,k) + Ts*eTrack;
    end

    if plotLqiLinearTest
        tTest = (0:N_lqiTest-1)*Ts;

        figure;
        plot(tTest, xTest.');
        grid on;
        xlabel('Time [s]');
        ylabel('State deviation');
        legend('theta1','theta2','theta1 dot','theta2 dot');
        title('LQI linear deviation test');

        figure;
        stairs(tTest(1:end-1), uCommandDevTest.');
        grid on;
        xlabel('Time [s]');
        ylabel('Command input deviation');
        title('LQI command deviation');
    end
end

%% Load EKF tuning result
ekfResultFile = fullfile(projectRoot, ...
    'extended_kalman_filter', 'ekf_tuning_result.mat');

ekfTuning = load(ekfResultFile, ...
    'Q_ekf', 'R_ekf', 'P0_ekf', 'Ts_ekf');

Q_ekf = ekfTuning.Q_ekf;
R_ekf = ekfTuning.R_ekf;
P0_ekf = ekfTuning.P0_ekf;
Ts_ekf = ekfTuning.Ts_ekf;

if abs(Ts_ekf - Ts) > 1e-12
    warning('EKF tuning sample time Ts_ekf = %.16g differs from plant Ts = %.16g.', Ts_ekf, Ts);
end

% Additional parameters/input vector for the EKF block, if you use a single
% vector input for Ts and p.
ekfInputParameters = [Ts; p(:)];

%% Compact diagnostics
fprintf('\nLQI/EKF setup complete.\n');
fprintf('Ts: %.6g s\n', Ts);
fprintf('Reference output: theta1 only\n');
fprintf('theta1_0: %.6f rad\n', theta1_0);
fprintf('Qx_lqi diag: ');
disp(diag(Qx_lqi).');
fprintf('Qi_lqi: %.6g\n', Qi_lqi);
fprintf('R_lqi: %.6g\n', R_lqi);
fprintf('Reference residual: %.3e\n', rel_res_ref);
fprintf('Largest |LQI pole|: %.6f\n', max(abs(closedLoopPoles_lqi)));

fprintf('\nKx_lqi_command:\n');
disp(Kx_lqi_command);
fprintf('Ki_lqi_command:\n');
disp(Ki_lqi_command);
fprintf('G_ref_command:\n');
disp(G_ref_command);

fprintf('\nUse in Simulink:\n');
fprintf('  e = r_dev - C_track*x_dev\n');
fprintf('  xi(k+1) = xi(k) + Ts*e(k)\n');
fprintf('  u_command_dev = Kx_lqi_command*x_dev + Ki_lqi_command*xi + G_ref_command*r_dev\n');
fprintf('  u_command = u0_command + u_command_dev\n');
