function [dx, y] = final_whitebox_7par_bias_fixed_ic_twofile_model(t, x, u, ...
    p_J1, p_kappa, p_g1, p_g2, p_Ku, p_b1, p_tau0, varargin)
%FINAL_WHITEBOX_7PAR_BIAS_FIXED_IC_TWOFILE_MODEL
%
% 7-parameter compact white-box model.
%
% States:
%   x(1) = theta_1       absolute rod-1 angle [rad]
%   x(2) = theta_2       absolute rod-2 angle [rad]
%   x(3) = omega_1       rod-1 angular velocity [rad/s]
%   x(4) = omega_2       rod-2 angular velocity [rad/s]
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
%   p_tau0   constant input/torque bias at joint 1
%
% Passive relative joint friction is fixed to zero.

theta1 = x(1);
theta2 = x(2);
omega1 = x(3);
omega2 = x(4);

if isempty(u)
    u_in = 0;
else
    u_in = u(1);
end

p_J1_safe = max(p_J1, 1e-6);

delta = theta1 - theta2;

% Inertia matrix. Positive-definite when p_J1 > 0 and |p_kappa| < 1.
coupling = p_kappa * sqrt(p_J1_safe);

M11 = p_J1_safe;
M12 = coupling * cos(delta);
M21 = M12;
M22 = 1.0;

% Velocity coupling terms.
C1 =  coupling * sin(delta) * omega2^2;
C2 = -coupling * sin(delta) * omega1^2;

% Gravity.
G1 = p_g1 * sin(theta1);
G2 = p_g2 * sin(theta2);

% Actuation and damping.
tau_motor = p_Ku * u_in + p_tau0;
F1 = p_b1 * omega1;

% Passive relative-joint friction removed.
F2 = 0.0;

rhs1 = tau_motor - F1 - G1 - C1 + F2;
rhs2 =             - G2 - C2 - F2;

detM = M11 * M22 - M12 * M21;

if abs(detM) < 1e-8
    detM = sign(detM + eps) * 1e-8;
end

alpha1 = ( M22 * rhs1 - M12 * rhs2) / detM;
alpha2 = (-M21 * rhs1 + M11 * rhs2) / detM;

dx = zeros(4, 1);
dx(1) = omega1;
dx(2) = omega2;
dx(3) = alpha1;
dx(4) = alpha2;

y = zeros(2, 1);
y(1) = theta1;
y(2) = theta2 - theta1;

end
