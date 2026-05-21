function [dx, y] = teammate_reduced_pendulum_model_offset(t, x, u, ...
    p_J, p_g1, p_g2, p_b1, p_b2, p_u, p_0, ...
    l1, g_const, theta1_offset, theta2_offset, varargin)
%TEAMMATE_REDUCED_PENDULUM_MODEL_OFFSET
% Teammate's reduced report model, implemented in a numerically stable form.
%
% States are LOCAL coordinates used in the identification data:
%   x(1) = theta_1_local  [rad]
%   x(2) = theta_2_local  [rad], relative passive joint angle
%   x(3) = omega_1        [rad/s]
%   x(4) = omega_2        [rad/s]
%
% The physical trigonometric angles are restored as:
%   theta_1_phys = theta_1_local + theta1_offset
%   theta_2_phys = theta_2_local + theta2_offset
%
% Your teammate's relation is used:
%   p_c = (l1/g_const) * p_g2
%
% Instead of estimating p_alpha directly, this file estimates
%   p_J = p_alpha - p_c^2.
% Then internally:
%   p_alpha = p_J + p_c^2.
% This is the same model family, but it guarantees a positive-definite
% inertia matrix when p_J > 0. For reporting, recover
%   p_alpha = p_J + p_c^2.
%
% No Coulomb or Stribeck friction is included.

%#ok<*INUSD>

theta1_local = x(1);
theta2_local = x(2);
omega1       = x(3);
omega2       = x(4);

if isempty(u)
    u_in = 0;
else
    u_in = u(1);
end

% Restore physical angles for sin/cos terms.
theta1 = theta1_local + theta1_offset;
theta2 = theta2_local + theta2_offset;

% Coupling parameter determined by p_g2 and known l1, g.
p_c = (l1 / g_const) * p_g2;
p_J_safe = max(p_J, 1e-8);
p_alpha = p_J_safe + p_c^2;

c = cos(theta2);
s = sin(theta2);

% Normalized inertia matrix from the reduced parameterisation.
M11 = p_alpha + 1 + 2*p_c*c;
M12 = 1 + p_c*c;
M22 = 1;

rhs1 = p_u*u_in + p_0 ...
     + p_c*s*(2*omega1*omega2 + omega2^2) ...
     - p_b1*omega1 ...
     + p_g1*sin(theta1) ...
     + p_g2*sin(theta1 + theta2);

rhs2 = -p_c*s*omega1^2 ...
     - p_b2*omega2 ...
     + p_g2*sin(theta1 + theta2);

detM = M11*M22 - M12*M12;
if abs(detM) < 1e-10
    detM = sign(detM + eps) * 1e-10;
end

theta1_ddot = ( M22*rhs1 - M12*rhs2) / detM;
theta2_ddot = (-M12*rhs1 + M11*rhs2) / detM;

dx = zeros(4, 1);
dx(1) = omega1;
dx(2) = omega2;
dx(3) = theta1_ddot;
dx(4) = theta2_ddot;

% Outputs are local measured coordinates.
y = zeros(2, 1);
y(1) = theta1_local;
y(2) = theta2_local;

end
