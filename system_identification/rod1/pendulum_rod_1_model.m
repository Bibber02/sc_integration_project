function [dx, y] = pendulum_rod_1_model( ...
    t, x, u, ...
    p_alpha, p_g1, p_b1, p_c1, eps_v1, p_u, p_0, l_1, ...
    p_b2, p_g2, p_c2, p_b_low, v_b, eps_v2, g, ...
    varargin)
%PENDULUM_ROD_1_MODEL
% Full nonlinear two-link rotational-pendulum model for idnlgrey, used
% to identify the rod-1-related (and motor-related) reduced parameters.
%
% Convention (same as the modelling chapter of the report):
%   theta_1 = 0  -> rod 1 pointing upward
%   theta_2 = 0  -> rod 2 aligned with the outward direction of rod 1
%   (theta_1, theta_2) = (pi, 0)  -> both rods hanging downward (stable)
%
% States:
%   x(1) = theta_1       [rad]
%   x(2) = theta_2       [rad]
%   x(3) = theta_1_dot   [rad/s]
%   x(4) = theta_2_dot   [rad/s]
%
% Input:
%   u    = motor command, scaled to [-1, +1]                  [-]
%
% Outputs:
%   y(1) = theta_1       [rad]
%   y(2) = theta_2       [rad]
%
% Parameters to estimate (rod 1 and motor; reduced/normalised by beta):
%   p_alpha  = alpha/beta                                      [-]
%   p_g1     = delta_1/beta = (m_1*c_1 + m_2*l_1)*g/beta       [1/s^2]
%   p_b1     = b_1_star/beta  (includes back-EMF)              [1/s]
%   p_c1     = F_C1/beta      (smooth Coulomb at joint 1)      [rad/s^2]
%   eps_v1   = smoothing velocity for joint-1 Coulomb          [rad/s]
%   p_u      = k_u/beta       (motor-command gain)             [1/s^2]
%   p_0      = tau_0/beta     (constant torque offset)         [1/s^2]
%   l_1      = rod-1 length (needed for p_c = (l_1/g)*p_g2)    [m]
%
% Auxiliary parameters (fixed; taken from rod-2 identification):
%   p_b2     = b_2/beta                                        [1/s]
%   p_g2     = delta_2/beta                                    [1/s^2]
%   p_c2     = F_C2/beta                                       [rad/s^2]
%   p_b_low  = low-velocity extra damping at joint 2           [1/s]
%   v_b      = low-velocity damping width                      [rad/s]
%   eps_v2   = smoothing velocity for joint-2 Coulomb          [rad/s]
%   g        = gravitational acceleration                      [m/s^2]
%
% The reduced equations of motion (see Section 4.2 of the report) are
%
%   (p_alpha + 1 + 2*p_c*cos(theta_2)) * theta_1_ddot
%   + (1 + p_c*cos(theta_2))          * theta_2_ddot
%   - p_c*sin(theta_2) * (2*theta_1_dot*theta_2_dot + theta_2_dot^2)
%   + F_1(theta_1_dot)/beta
%   - p_g1*sin(theta_1)
%   - p_g2*sin(theta_1 + theta_2)
%   = p_u*u + p_0,
%
%   (1 + p_c*cos(theta_2)) * theta_1_ddot
%   +                       theta_2_ddot
%   + p_c*sin(theta_2)*theta_1_dot^2
%   + F_2(theta_2_dot)/beta
%   - p_g2*sin(theta_1 + theta_2)
%   = 0,
%
% with smooth (already-divided-by-beta) friction laws
%
%   F_1(theta_1_dot)/beta =  p_b1*theta_1_dot
%                          + p_c1*tanh(theta_1_dot/eps_v1),
%
%   F_2(theta_2_dot)/beta =  [p_b2 + p_b_low*exp(-(theta_2_dot/v_b)^2)]
%                            * theta_2_dot
%                          + p_c2*tanh(theta_2_dot/eps_v2).
%
% The geometric coupling constant p_c is computed from p_g2 and l_1:
%
%   p_c = (l_1/g) * p_g2.

% ------------------------------------------------------------------
% Unpack states
% ------------------------------------------------------------------
theta_1     = x(1);
theta_2     = x(2);
theta_1_dot = x(3);
theta_2_dot = x(4);

% ------------------------------------------------------------------
% Guard the smoothing parameters against zero (numerical safety)
% ------------------------------------------------------------------
eps_v1 = max(eps_v1, 1e-8);
eps_v2 = max(eps_v2, 1e-8);
v_b    = max(v_b,    1e-8);

% ------------------------------------------------------------------
% Derived geometric coupling: p_c = (l_1/g) * p_g2
% ------------------------------------------------------------------
p_c = (l_1 / g) * p_g2;

% ------------------------------------------------------------------
% Friction torques (already normalised by beta)
% ------------------------------------------------------------------
F1_over_beta =  p_b1 * theta_1_dot ...
              + p_c1 * tanh(theta_1_dot / eps_v1);

low_velocity_multiplier = exp(-(theta_2_dot / v_b)^2);
effective_b2 = p_b2 + p_b_low * low_velocity_multiplier;

F2_over_beta =  effective_b2 * theta_2_dot ...
              + p_c2 * tanh(theta_2_dot / eps_v2);

% ------------------------------------------------------------------
% Pre-compute trig quantities
% ------------------------------------------------------------------
c2 = cos(theta_2);
s2 = sin(theta_2);
s1   = sin(theta_1);
s12  = sin(theta_1 + theta_2);

% ------------------------------------------------------------------
% Mass matrix M_bar = M / beta  (symmetric 2x2)
% ------------------------------------------------------------------
M11 = p_alpha + 1 + 2*p_c*c2;
M12 = 1 + p_c*c2;
M22 = 1;

% ------------------------------------------------------------------
% Right-hand side b_bar = (tau_ext - C*qdot - F - g_vec) / beta
% Built from the two scalar EoMs by moving everything except the
% q_ddot terms to the right-hand side.
% ------------------------------------------------------------------
rhs1 =   p_u*u + p_0 ...
       + p_c*s2 * (2*theta_1_dot*theta_2_dot + theta_2_dot^2) ...
       - F1_over_beta ...
       + p_g1*s1 ...
       + p_g2*s12;

rhs2 = - p_c*s2 * theta_1_dot^2 ...
       - F2_over_beta ...
       + p_g2*s12;

% ------------------------------------------------------------------
% Solve M_bar * [theta_1_ddot; theta_2_ddot] = [rhs1; rhs2]
% explicitly (2x2 inverse; det = M11*M22 - M12^2 = M11 - M12^2).
% ------------------------------------------------------------------
detM = M11*M22 - M12*M12;

theta_1_ddot = ( M22*rhs1 - M12*rhs2 ) / detM;
theta_2_ddot = (-M12*rhs1 + M11*rhs2 ) / detM;

% ------------------------------------------------------------------
% Pack state derivative
% ------------------------------------------------------------------
dx = zeros(4,1);
dx(1) = theta_1_dot;
dx(2) = theta_2_dot;
dx(3) = theta_1_ddot;
dx(4) = theta_2_ddot;

% ------------------------------------------------------------------
% Outputs: both angles are measured
% ------------------------------------------------------------------
y = [theta_1; theta_2];

end