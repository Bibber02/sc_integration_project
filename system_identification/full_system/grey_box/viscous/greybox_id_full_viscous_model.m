function [dx, y] = greybox_id_full_viscous_model(t, x, u, p_a, p_b1, p_g1, p_u, p_0, p_b2, p_g2, varargin)
%GREYBOX_ID_FULL_VISCOUS_MODEL Reduced full-system grey-box model.
%
% States:
%   x(1) = theta_1       rod-1 angle [rad]
%   x(2) = theta_2       rod-2 relative angle [rad]
%   x(3) = theta_1_dot   rod-1 angular velocity [rad/s]
%   x(4) = theta_2_dot   rod-2 relative angular velocity [rad/s]
%
% Outputs:
%   y(1) = theta_1
%   y(2) = theta_2
%
% Friction:
%   Joint 1: viscous only
%   Joint 2: viscous only

% Known constants
g  = 9.81;   % [m/s^2]
l1 = 0.10;   % [m]

% Dependent reduced coupling parameter
p_c = l1 * p_g2 / g;

% States
theta_1     = x(1);
theta_2     = x(2);
theta_1_dot = x(3);
theta_2_dot = x(4);

% Input voltage
u_1 = u(1);

% Trigonometric terms
cos_theta_2 = cos(theta_2);
sin_theta_2 = sin(theta_2);

% Viscous friction only
F_1 = p_b1 * theta_1_dot;
F_2 = p_b2 * theta_2_dot;

% Reduced inertia matrix M_r(q)
M_11 = p_a + 1 + 2 * p_c * cos_theta_2;
M_12 = 1 + p_c * cos_theta_2;
M_21 = M_12;
M_22 = 1;

M = [M_11, M_12;
     M_21, M_22];

% Reduced Coriolis/centrifugal vector C_r(q,qdot) qdot
C_11 = -p_c * sin_theta_2 * theta_2_dot;
C_12 = -p_c * sin_theta_2 * (theta_1_dot + theta_2_dot);
C_21 =  p_c * sin_theta_2 * theta_1_dot;
C_22 = 0;

C_vec = [C_11 * theta_1_dot + C_12 * theta_2_dot;
         C_21 * theta_1_dot + C_22 * theta_2_dot];

% Reduced gravity vector
G = [-p_g1 * sin(theta_1) - p_g2 * sin(theta_1 + theta_2);
     -p_g2 * sin(theta_1 + theta_2)];

% Input vector: motor voltage gain plus constant bias
I = [p_u * u_1 + p_0;
     0];

% Dynamics:
%   M*qddot + C*qdot + F + G = I
rhs = I - C_vec - [F_1; F_2] - G;
theta_ddot = M \ rhs;

% State derivative
dx = zeros(4, 1);
dx(1) = theta_1_dot;
dx(2) = theta_2_dot;
dx(3) = theta_ddot(1);
dx(4) = theta_ddot(2);

% Output vector
y = zeros(2, 1);
y(1) = theta_1;
y(2) = theta_2;

end
