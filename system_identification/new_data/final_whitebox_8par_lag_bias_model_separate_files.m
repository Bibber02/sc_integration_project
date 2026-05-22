function [dx, y] = final_whitebox_8par_lag_bias_model_separate_files(t, x, u, ...
    p_J1, p_kappa, p_g1, p_g2, p_Ku, p_b1, p_tau0, p_Tm, varargin)
%FINAL_WHITEBOX_8PAR_LAG_BIAS_MODEL_SEPARATE_FILES
%
% Compact white-box / grey-box model for the whole system.
%
% This version is meant for separate training and validation recordings.
%
% Model choices:
%   - absolute angles internally;
%   - relative second encoder output;
%   - no passive rod-2 friction term;
%   - effective viscous damping only at the actuated first joint;
%   - first-order actuator lag;
%   - constant input/torque bias term.
%
% States:
%   x(1) = theta_1       absolute rod-1 angle [rad]
%   x(2) = theta_2       absolute rod-2 angle [rad]
%   x(3) = omega_1       rod-1 angular velocity [rad/s]
%   x(4) = omega_2       rod-2 angular velocity [rad/s]
%   x(5) = u_act         filtered/effective actuator input [-]
%
% Outputs:
%   y(1) = theta_1
%   y(2) = phi_2 = theta_2 - theta_1
%
% Estimated parameters:
%   p_J1     rod-1/motor-side normalized inertia
%   p_kappa  normalized coupling coefficient, constrained to |p_kappa| < 1
%   p_g1     rod-1 gravity coefficient
%   p_g2     rod-2 gravity coefficient
%   p_Ku     input-to-torque gain
%   p_b1     effective rod-1 damping, including motor back-EMF
%   p_tau0   constant torque/input bias at joint 1
%   p_Tm     first-order actuator lag time constant [s]
%
% Passive relative joint friction is fixed to zero:
%   F2 = 0

% States
theta1 = x(1);
theta2 = x(2);
omega1 = x(3);
omega2 = x(4);
u_act  = x(5);

if isempty(u)
    u_in = 0;
else
    u_in = u(1);
end

% Numerical guards
p_J1_safe = max(p_J1, 1e-6);
p_Tm_safe = max(p_Tm, 1e-4);

delta = theta1 - theta2;

% Inertia matrix.
% Off-diagonal term is parameterized as p_kappa*sqrt(p_J1), which helps keep
% the inertia matrix positive definite when p_J1 > 0 and |p_kappa| < 1.
coupling = p_kappa * sqrt(p_J1_safe);

M11 = p_J1_safe;
M12 = coupling * cos(delta);
M21 = M12;
M22 = 1.0;

% Velocity coupling terms.
C1 =  coupling * sin(delta) * omega2^2;
C2 = -coupling * sin(delta) * omega1^2;

% Gravity terms.
G1 = p_g1 * sin(theta1);
G2 = p_g2 * sin(theta2);

% Actuated joint damping and motor input.
F1 = p_b1 * omega1;
tau_motor = p_Ku * u_act + p_tau0;

% Passive relative-joint damping/friction removed in this compact model.
F2 = 0.0;

% Generalized torque balance.
rhs1 = tau_motor - F1 - G1 - C1 + F2;
rhs2 =             - G2 - C2 - F2;

% Solve M*qddot = rhs.
detM = M11 * M22 - M12 * M21;

if abs(detM) < 1e-8
    detM = sign(detM + eps) * 1e-8;
end

alpha1 = ( M22 * rhs1 - M12 * rhs2) / detM;
alpha2 = (-M21 * rhs1 + M11 * rhs2) / detM;

% First-order actuator lag.
u_act_dot = (u_in - u_act) / p_Tm_safe;

% State derivatives
dx = zeros(5, 1);
dx(1) = omega1;
dx(2) = omega2;
dx(3) = alpha1;
dx(4) = alpha2;
dx(5) = u_act_dot;

% Relative encoder output for rod 2.
y = zeros(2, 1);
y(1) = theta1;
y(2) = theta2 - theta1;

end
