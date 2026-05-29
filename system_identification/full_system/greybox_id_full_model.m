function [dx, y] = greybox_id_full_model(t, x, u, p_a, p_b1, p_c1, p_g1, p_u, p_0, p_b2, p_g2, p_c2, p_sdelta2, v_s2, eps_v1, eps_v2, varargin)

% States:
%   x(1) = theta_1       relative rod-1 angle [rad]
%   x(2) = theta_2       relative rod-2 angle [rad]
%   x(3) = omega_1       rod-1 angular velocity [rad/s]
%   x(4) = omega_2       rod-2 angular velocity [rad/s]
%
% Outputs:
%   y(1) = theta_1
%   y(2) = theta_2

% Known constants
g  = 9.81;   % [m/s^2]
l1 = 0.10;   % [m]

% Dependent coupling parameter
p_c = l1 * p_g2 / g;

% States
theta1 = x(1);
theta2 = x(2);
theta1_dot = x(3);
theta2_dot = x(4);

% Input
u1 = u(1);

% Precompute trigonometric terms
c2 = cos(theta2);
s2 = sin(theta2);

% --- Friction Models ---
% Active joint 1 (Coulomb + Viscous)
eps_v1 = max(eps_v1, 1e-5);
F1 = p_b1 * theta1_dot + p_c1 * tanh(theta1_dot / eps_v1);

% Passive joint 2 (Stribeck from report)
v_s2 = max(v_s2, 1e-5);
eps_v2 = max(eps_v2, 1e-5);
friction_level2 = p_c2 + p_sdelta2 * exp(-(theta2_dot / v_s2)^2);
F2 = p_b2 * theta2_dot + friction_level2 * tanh(theta2_dot / eps_v2);

% --- Dynamics Matrices ---
% Inertia matrix M_r(q)
M11 = p_a + 1 + 2*p_c*c2;
M12 = 1 + p_c*c2;
M21 = M12;
M22 = 1.0;
M = [M11 M12; 
     M21 M22];

% Coriolis/centrifugal matrix C_r(q,qdot) * qdot
C11 = -p_c * s2 * theta2_dot;
C12 = -p_c * s2 * (theta1_dot + theta2_dot);
C21 =  p_c * s2 * theta1_dot;
C22 = 0;
C_vec = [C11 * theta1_dot + C12 * theta2_dot;
         C21 * theta1_dot + C22 * theta2_dot];

% Gravity vector g_r(q)
G = [-p_g1 * sin(theta1) - p_g2 * sin(theta1 + theta2);
     -p_g2 * sin(theta1 + theta2)];

% Input vector (Motor + Torque Offset)
I = [p_u * u1 + p_0;
     0];

% Solve: M*qddot + C_vec + F + G = I
rhs = I - C_vec - [F1; F2] - G;
theta_dd = M \ rhs;

% State derivative
dx = zeros(4, 1);
dx(1) = theta1_dot;
dx(2) = theta2_dot;
dx(3) = theta_dd(1);
dx(4) = theta_dd(2);

% Output vector
y = zeros(2, 1);
y(1) = theta1;
y(2) = theta2;

end