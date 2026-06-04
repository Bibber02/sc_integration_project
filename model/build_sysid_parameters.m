%% setup_pendulum_model.m
% Loads identified parameters, physical constants, derived quantities and
% default initial conditions into the base workspace so that the Simulink
% plant model of the rotational pendulum can pick them up by name.
%
% Convention (see report, Section 3.1):
%   theta_1 : angle of rod 1, measured from the UPWARD vertical (CW +)
%   theta_2 : angle of rod 2, measured from the OUTWARD direction of rod 1
%   (theta_1, theta_2) = (0,0)   -> fully inverted (both rods up)
%   (theta_1, theta_2) = (pi,0)  -> stable hanging equilibrium
%
% Model form (reduced/normalised by beta = I_2 + m_2 c_2^2):
%   M(q) qddot + C(q,qdot) qdot + F(qdot) + G(q) = I(u)
%
% All parameters below are the "reduced" parameters used during the
% grey-box identification.

clear; clc;

%% ----- Identified parameters --------------------------------
% Rod 1
p_alpha = 33.6175;        % [-]       Normalised inertia ratio alpha/beta
p_b1    = 744.5341;         % [1/s]     Viscous damping at joint 1 (incl. back-EMF)
p_c1    = 23.9089;        % [rad/s^2] Kinetic Coulomb friction strength
p_g1    = 105.6867;        % [1/s^2]   Normalised gravity parameter at joint 1
p_u     = 4091.4423;        % [1/s^2]   motor-command gain ku/beta
p_0     = -40.8063;       % [1/s^2]   constant torque offset tau_0/beta
p_0 = 0;

% Rod 2
p_b2    = 0.0670;    % [1/s]     viscous friction at joint 2;
p_g2    = 112.0267;        % [1/s^2]   normalised gravity coefficient at joint 2
p_c2    = 0.2491;       % [rad/s^2] Kinetic Coulomb friction at joint 2
p_s_d2  = 0.1836;       % [1/s]     Static-minus-Coulomb friction strength at Join 2
v_s2    = 3.5000;       % [rad/s]   Stribeck Transition Velocity

%% ----- Fixed (non-identified) smoothing & physical constants
eps_v1 = 0.02;            % [rad/s]   tanh-smoothing velocity, joint 1
eps_v2 = 0.02;            % [rad/s]   tanh-smoothing velocity, joint 2
l_1    = 0.100;           % [m]       motor pivot to swivel joint
g      = 9.81;            % [m/s^2]   gravitational acceleration

%% ----- Derived constants ------------------------------------
% Geometric coupling parameter
p_c            = (l_1 / g) * p_g2;   % [-]

%% ----- Simulation / sampling settings -----------------------
h      = 0.01;            % [s]       sampling period (matches real-time)

%% ----- Default initial conditions ---------------------------
% Stable hanging equilibrium by default. Switch to (0,0) for swing-up.
theta_1_0     = pi;       % [rad]
theta_2_0     = 0;        % [rad]
theta_1_dot_0 = 0;        % [rad/s]
theta_2_dot_0 = 0;        % [rad/s]

x0 = [theta_1_0; theta_2_0; theta_1_dot_0; theta_2_dot_0];

%% ----- Sanity print -----------------------------------------
fprintf('Rotational pendulum parameters loaded.\n');
fprintf('  Reduced inertia coupling p_c = %.4f\n', p_c);
fprintf('  Sampling period h            = %.3f s\n', h);
fprintf('  Initial state x0             = [%.3f %.3f %.3f %.3f]^T\n', x0);