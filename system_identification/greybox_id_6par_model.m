function [dx, y] = greybox_id_6par_model(t, x, u, p_a, p_b1, p_b2, p_g1, p_g2, p_u, varargin)

% States:
%   x(1) = theta_1       relative rod-1 angle [rad]
%   x(2) = theta_2       relative rod-2 angle [rad]
%   x(3) = omega_1       rod-1 angular velocity [rad/s]
%   x(4) = omega_2       rod-2 angular velocity [rad/s]
%
% Outputs:
%   y(1) = theta_1
%   y(2) = theta_2
%
% Estimated parameters:
%   p_a
%   p_b1
%   p_b2
%   p_g1
%   p_g2
%   p_u

% Known constants
g  = 9.81;   % [m/s^2]
l1 = 0.10;   % [m]

% Dependent coupling parameter
p_c = l1 * p_g2 / g;

% States
theta1     = x(1);
theta2     = x(2);
theta1_dot = x(3);
theta2_dot = x(4);

% Input
u1 = u(1);

% Precompute trigonometric terms
c2 = cos(theta2);
s2 = sin(theta2);

% Inertia matrix M_r(q)
M11 = p_a + 1 + 2*p_c*c2;
M12 = 1 + p_c*c2;
M21 = M12;
M22 = 1.0;

M = [M11 M12;
     M21 M22];

% Coriolis/centrifugal matrix C_r(q,qdot)
C11 = -p_c * s2 * theta2_dot;
C12 = -p_c * s2 * (theta1_dot + theta2_dot);
C21 =  p_c * s2 * theta1_dot;
C22 = 0;

C = [C11 C12;
     C21 C22];

% Viscous damping matrix D_r
D = [p_b1 0;
     0    p_b2];

% Gravity vector g_r(q)
G = [-p_g1 * sin(theta1) - p_g2 * sin(theta1 + theta2);
     -p_g2 * sin(theta1 + theta2)];

% Input vector b_r u
I = [p_u * u1;
     0];

% Generalized velocity vector
theta_d = [theta1_dot;
           theta2_dot];

% Solve:
% M*qddot + C*qdot + D*qdot + G = I
rhs = I - C*theta_d - D*theta_d - G;

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
