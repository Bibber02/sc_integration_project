clear;
clc;

%% User settings
sampleTime = 0.01;

% Equilibrium:
%   [0;  0; 0; 0] = both rods up
%   [pi; 0; 0; 0] = both rods down
%   [pi; pi; 0; 0] = rod 1 down, rod 2 up
x0 = [pi; pi; 0; 0];

Q_lqr = diag([10 10 0.1 0.1]);
R_lqr = 1;

%% Kalman filter settings
% These matrices are used by the Simulink Kalman Filter block.
Q_kf = diag([5e-4 3e-4 1e-6 1e-6])*0.08;
R_kf = diag([1e-4 1e-4]);

% Initial covariance.
P0_kf = diag([1e-3 1e-3 1e-1 1e-1]);

% The linearized model is a deviation model, so the KF state also starts
% from deviation coordinates.
xhat0_kf = zeros(4,1);

%% Paths
scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptFolder));

modelFolder = fullfile(projectRoot, 'model');
greyBoxFolder = fullfile(modelFolder, 'grey_box_stribeck');
hardwareFolder = fullfile(scriptFolder, 'rotating-pendulum');

addpath(modelFolder, '-begin');
addpath(greyBoxFolder, '-begin');
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

% Useful constants for Simulink deviation signals
x0_kf = x0;
u0_kf = u0;
y0_kf = y0;

[K_lqr, closedLoopPoles] = calc_lqr(sys_disc, Q_lqr, R_lqr);

%% LQR poles
s_closedLoopPoles = log(closedLoopPoles)/Ts;

fprintf('\nLQR closed-loop poles in z-domain:\n');
disp(closedLoopPoles);

fprintf('LQR closed-loop poles in equivalent s-domain:\n');
disp(s_closedLoopPoles);

%% Kalman filter matrices for Simulink
A_kf = Ad;
B_kf = Bd;
C_kf = Cd;
D_kf = Dd;

% Process noise enters all states directly
G_kf = eye(size(A_kf,1));

% Steady-state Kalman gain
[L_kf, P_kf] = dlqe(A_kf, G_kf, C_kf, Q_kf, R_kf);

% Correct observer pole calculation
observerPoles = eig(A_kf - L_kf*C_kf);

% Equivalent continuous-time poles
s_observerPoles = log(observerPoles)/Ts;

fprintf('\nKalman observer poles in z-domain:\n');
disp(observerPoles);

fprintf('Kalman observer poles in equivalent s-domain:\n');
disp(s_observerPoles);

fprintf('Largest observer pole magnitude: %.6f\n', max(abs(observerPoles)));

disp('Kalman gain L_kf = ');
disp(L_kf);

%% Reference tracking gain for theta1 only
% Track:
%
%   z = theta1
%
% Controller:
%
%   u_dev = -K_lqr*x_hat_dev + G_ref*r_dev

Ctilde = [1 0 0 0];
Dtilde = 0;

n = size(Ad,1);
m = size(Bd,2);

M_ref = [Ad - eye(n), Bd;
         Ctilde,      Dtilde];

rhs_ref = [zeros(n,1);
           1];

X_ref = M_ref \ rhs_ref;

Pi_ref    = X_ref(1:n,:);
Gamma_ref = X_ref(n+1:n+m,:);

G_ref = Gamma_ref + K_lqr*Pi_ref;

% Equilibrium theta1 around which the linear model was made
theta1_0 = Ctilde*x0 + Dtilde*u0;

% Check residual
res_ref = M_ref*X_ref - rhs_ref;
rel_res_ref = norm(res_ref) / max(1,norm(rhs_ref));

fprintf('\nReference tracking gain for theta1:\n');
fprintf('G_ref = %.6f\n', G_ref);
fprintf('theta1_0 = %.6f rad\n', theta1_0);
fprintf('Relative regulator-equation residual: %.3e\n', rel_res_ref);

fs = 1/h;
fc = 10;
wn = fc / (fs/2);
[b, a] = butter(1, wn);