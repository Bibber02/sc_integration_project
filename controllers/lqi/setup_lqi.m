%% Simple LQI + EKF setup for selected no-p0 down-line model
clear;
clc;

%% User settings
sampleTime = 0.01;

% Sign convention between model input and actual command input.
%   u_model_dev = inputSignCommandToModel * u_command_dev
% Use -1 if your physical command sign is opposite to the identified model.
inputSignCommandToModel = -1;

% Determines what output we want the integral action to track.
% Current choice: theta1 only, in measured/saved coordinates.
C_track = [1 0 0 0];
D_track = 0;

% a_lqr penalizes active angle deviation.
% b_lqr penalizes passive absolute angle deviation.
%
% Since theta2 is relative, the passive absolute angle deviation is
% approximately theta1_dev + theta2_dev in the measured-coordinate deviation
% model. The common theta_scale only multiplies this cost and therefore does
% not change the structure of the Q block.
a_lqr = 5;
b_lqr = 100;

Qx_lqi = [a_lqr+b_lqr, b_lqr, 0, 0;
          b_lqr,       b_lqr, 0, 0;
          0,           0,     0, 0;
          0,           0,     0, 0];

% These parameters control the cost of integral error and control effort.
Qi_lqi = 5;
R_lqi  = 1;

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

ekfStateFcn = str2func(ekfStateFcnName);
ekfMeasurementFcn = str2func(ekfMeasurementFcnName);

%% Hardware initialization
run(fullfile(hardwareFolder, 'hwinit.m'));

%% Plant model
[p, parameterInfo] = load_parameters();

% Use the measured-coordinate state corresponding to:
%   theta1_phys = pi, theta2_phys = 0, velocities = 0.
% This is not necessarily exactly [pi; 0; 0; 0] in measured coordinates
% once the down-line calibration correction is included.
x0 = measuredAllDownEquilibriumFromParameters(p);

lin = linearize_rotpendulum(struct( ...
    'x0', x0, ...
    'u0', 0, ...
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
% These are useful in Simulink for subtracting the operating point and
% adding the equilibrium input back to the controller output.
x_eq_lqi = x0;
y_eq_lqi = y0;
u0_model = u0;
u0_command = inputSignCommandToModel * u0_model;

% Initial/reference offset for the tracked output.
theta1_0 = C_track*x_eq_lqi + D_track*u0_model;
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

Acl_aug_lqi = A_aug - B_aug*K_aug;

closedLoopPoles_lqi = eig(Acl_aug_lqi);
s_closedLoopPoles_lqi = log(closedLoopPoles_lqi)/Ts;

%% Closed-loop damping diagnostics
z_poles_lqi = closedLoopPoles_lqi;
s_poles_lqi = s_closedLoopPoles_lqi;

zeta_lqi = nan(size(s_poles_lqi));
omega_n_lqi = nan(size(s_poles_lqi));

zPoleImagTol_lqi = 1e-8;
isComplexMode_lqi = abs(imag(z_poles_lqi)) > zPoleImagTol_lqi;
isNegativeRealZ_lqi = ~isComplexMode_lqi & real(z_poles_lqi) < 0;

omega_n_lqi(isComplexMode_lqi) = abs(s_poles_lqi(isComplexMode_lqi));
zeta_lqi(isComplexMode_lqi) = -real(s_poles_lqi(isComplexMode_lqi)) ./ ...
                               omega_n_lqi(isComplexMode_lqi);

decayRate_lqi = log(abs(z_poles_lqi))/Ts;

[~, dominantPoleOrder_lqi] = sort(abs(z_poles_lqi), 'descend');

complexPoleOrder_lqi = dominantPoleOrder_lqi(isComplexMode_lqi(dominantPoleOrder_lqi));
if ~isempty(complexPoleOrder_lqi)
    dominantComplexPoleIndex_lqi = complexPoleOrder_lqi(1);
    dominantComplexPole_lqi = z_poles_lqi(dominantComplexPoleIndex_lqi);
    dominantComplexPole_s_lqi = s_poles_lqi(dominantComplexPoleIndex_lqi);
    dominantComplexZeta_lqi = zeta_lqi(dominantComplexPoleIndex_lqi);
    dominantComplexOmegaN_lqi = omega_n_lqi(dominantComplexPoleIndex_lqi);
else
    dominantComplexPoleIndex_lqi = [];
    dominantComplexPole_lqi = [];
    dominantComplexPole_s_lqi = [];
    dominantComplexZeta_lqi = NaN;
    dominantComplexOmegaN_lqi = NaN;
end

%% Reference feedforward matrix
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
Kx_lqi_model = Kx_lqi;
Ki_lqi_model = Ki_lqi;

Kx_lqi_command = -inputSignCommandToModel * Kx_lqi_model;
Ki_lqi_command = -inputSignCommandToModel * Ki_lqi_model;

K_lqi_model = [Kx_lqi_model Ki_lqi_model];
K_lqi_command = [Kx_lqi_command Ki_lqi_command];

%% Load EKF tuning result
ekfResultFile = fullfile(projectRoot, ...
    'extended_kalman_filter', 'ekf_tuning_result.mat');

ekfTuning = load(ekfResultFile, ...
    'Q_ekf', 'R_ekf', 'P0_ekf', 'Ts_ekf');

Q_ekf = ekfTuning.Q_ekf;
R_ekf = ekfTuning.R_ekf;
P0_ekf = ekfTuning.P0_ekf;
Ts_ekf = ekfTuning.Ts_ekf;

% Additional parameters/input vector for the EKF block when using one muxed
% input. Feed the EKF transition function [u; Ts; p(:)].
ekfInputParameters = [Ts; p(:)];

%% Compact diagnostics
thetaScale = p(13);
theta1Offset = p(14);
thetaAbsDownRaw = p(15);
theta2Offset = pi - thetaScale*thetaAbsDownRaw - theta1Offset;

theta1PhysEq = thetaScale*x_eq_lqi(1) + theta1Offset;
theta2PhysEq = thetaScale*x_eq_lqi(2) + theta2Offset;

fprintf('\nLQI/EKF setup complete.\n');
fprintf('Selected parameter file:\n  %s\n', parameterInfo.matFile);
fprintf('Ts: %.6g s\n', Ts);
fprintf('Reference output: theta1 only\n');
fprintf('Measured-coordinate equilibrium x_eq_lqi:\n');
disp(x_eq_lqi);
fprintf('Physical equilibrium check: theta1_phys = %.9g rad, theta2_phys = %.9g rad, sum = %.9g rad\n', ...
    theta1PhysEq, theta2PhysEq, theta1PhysEq + theta2PhysEq);
fprintf('theta1_0: %.6f rad\n', theta1_0);
fprintf('u0_model: %.6f\n', u0_model);
fprintf('u0_command: %.6f\n', u0_command);
fprintf('f0 at linearization point:\n');
disp(f0);
fprintf('a_lqr: %.6g\n', a_lqr);
fprintf('b_lqr: %.6g\n', b_lqr);
fprintf('Qx_lqi:\n');
disp(Qx_lqi);
fprintf('Qi_lqi: %.6g\n', Qi_lqi);
fprintf('R_lqi: %.6g\n', R_lqi);
fprintf('Reference residual: %.3e\n', rel_res_ref);
fprintf('Largest |LQI pole|: %.6f\n', max(abs(closedLoopPoles_lqi)));

fprintf('\nClosed-loop poles sorted by |z|:\n');
fprintf('  index        |z|          z pole                    s pole [1/s]              zeta / note        omega_n / decay\n');
for kk = 1:length(dominantPoleOrder_lqi)
    ii = dominantPoleOrder_lqi(kk);

    if isComplexMode_lqi(ii)
        fprintf('  %2d      %.6f    %.6f%+.6fi    %.4f%+.4fi    zeta = %.4f      wn = %.4f\n', ...
            ii, abs(z_poles_lqi(ii)), ...
            real(z_poles_lqi(ii)), imag(z_poles_lqi(ii)), ...
            real(s_poles_lqi(ii)), imag(s_poles_lqi(ii)), ...
            zeta_lqi(ii), omega_n_lqi(ii));

    elseif isNegativeRealZ_lqi(ii)
        fprintf('  %2d      %.6f    %.6f%+.6fi    %.4f%+.4fi    negative real z   decay = %.4f\n', ...
            ii, abs(z_poles_lqi(ii)), ...
            real(z_poles_lqi(ii)), imag(z_poles_lqi(ii)), ...
            real(s_poles_lqi(ii)), imag(s_poles_lqi(ii)), ...
            decayRate_lqi(ii));

    else
        fprintf('  %2d      %.6f    %.6f%+.6fi    %.4f%+.4fi    real mode         decay = %.4f\n', ...
            ii, abs(z_poles_lqi(ii)), ...
            real(z_poles_lqi(ii)), imag(z_poles_lqi(ii)), ...
            real(s_poles_lqi(ii)), imag(s_poles_lqi(ii)), ...
            decayRate_lqi(ii));
    end
end

if ~isempty(dominantComplexPoleIndex_lqi)
    fprintf('\nDominant complex mode: pole index %d, zeta = %.4f, omega_n = %.4f rad/s\n', ...
        dominantComplexPoleIndex_lqi, dominantComplexZeta_lqi, dominantComplexOmegaN_lqi);
else
    fprintf('\nNo complex closed-loop poles found, so no oscillatory damping ratio was calculated.\n');
end

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

%% Local helper
function x_eq = measuredAllDownEquilibriumFromParameters(p)
theta_scale        = max(p(13), 1e-6);
theta1_offset      = p(14);
theta_abs_down_raw = p(15);
theta2_offset      = pi - theta_scale*theta_abs_down_raw - theta1_offset;

theta1_meas_eq = (pi - theta1_offset) / theta_scale;
theta2_meas_eq = (0  - theta2_offset) / theta_scale;

x_eq = [theta1_meas_eq; theta2_meas_eq; 0; 0];
end
