function [xdot, y] = pendulum_nonlinear_rhs(x, u, p)
%PENDULUM_NONLINEAR_RHS Nonlinear identified plant model for Simulink.
%
% States:
%   x(1) = theta_1
%   x(2) = theta_2
%   x(3) = theta_1_dot
%   x(4) = theta_2_dot
%
% Input:
%   u = motor voltage / command
%
% Output:
%   y = [theta_1; theta_2]

% Unpack parameters
p_a        = p(1);
p_b1       = p(2);
p_c1       = p(3);
p_g1       = p(4);
p_u        = p(5);
p_0        = p(6);
p_b2       = p(7);
p_g2       = p(8);
p_c2       = p(9);
p_sdelta2  = p(10);
v_s2       = p(11);
eps_v1     = p(12);
eps_v2     = p(13);

% Known constants
g  = 9.81;
l1 = 0.10;

% Dependent coupling parameter
p_c = l1 * p_g2 / g;

% States
theta1     = x(1);
theta2     = x(2);
theta1_dot = x(3);
theta2_dot = x(4);

% Input
u1 = u(1);

% Trigonometry
c2 = cos(theta2);
s2 = sin(theta2);

% Friction joint 1
eps_v1 = max(eps_v1, 1e-5);
F1 = p_b1 * theta1_dot + p_c1 * tanh(theta1_dot / eps_v1);

% Friction joint 2
v_s2 = max(v_s2, 1e-5);
eps_v2 = max(eps_v2, 1e-5);
friction_level2 = p_c2 + p_sdelta2 * exp(-(theta2_dot / v_s2)^2);
F2 = p_b2 * theta2_dot + friction_level2 * tanh(theta2_dot / eps_v2);

% Inertia matrix
M11 = p_a + 1 + 2*p_c*c2;
M12 = 1 + p_c*c2;
M21 = M12;
M22 = 1.0;

M = [M11 M12;
     M21 M22];

% Coriolis/centrifugal vector
C11 = -p_c * s2 * theta2_dot;
C12 = -p_c * s2 * (theta1_dot + theta2_dot);
C21 =  p_c * s2 * theta1_dot;
C22 = 0;

C_vec = [C11 * theta1_dot + C12 * theta2_dot;
         C21 * theta1_dot + C22 * theta2_dot];

% Gravity vector
G = [-p_g1 * sin(theta1) - p_g2 * sin(theta1 + theta2);
     -p_g2 * sin(theta1 + theta2)];

% Input vector
I = [p_u * u1;
     0];

% Dynamics
rhs = I - C_vec - [F1; F2] - G;
theta_dd = M \ rhs;

% State derivative
xdot = zeros(4,1);
xdot(1) = theta1_dot;
xdot(2) = theta2_dot;
xdot(3) = theta_dd(1);
xdot(4) = theta_dd(2);

% Output
y = zeros(2,1);
y(1) = theta1;
y(2) = theta2;

end