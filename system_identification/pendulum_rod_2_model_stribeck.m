function [dx, y] = pendulum_rod_2_model_stribeck( ...
    t, x, u, p_b2, p_g2, p_c2, p_delta_s2, v_s2, eps_v, varargin)
%PENDULUM_ROD_2_MODEL_STRIBECK Nonlinear passive-link model for idnlgrey.
%
% Model:
%
%   theta_ddot =
%       -p_b2*theta_dot
%       -[p_c2 + p_delta_s2*exp(-(theta_dot/v_s2)^2)]
%           *tanh(theta_dot/eps_v)
%       -p_g2*sin(theta)
%
% States:
%   x(1) = theta      [rad]
%   x(2) = theta_dot  [rad/s]
%
% Parameters:
%   p_b2        = viscous damping parameter          [1/s]
%   p_g2        = gravity parameter                  [1/s^2]
%   p_c2        = Coulomb friction strength          [rad/s^2]
%   p_delta_s2  = static-minus-Coulomb friction      [rad/s^2]
%   v_s2        = Stribeck transition velocity       [rad/s]
%   eps_v       = smoothing velocity                 [rad/s]
%
% Static friction strength:
%
%   p_s2 = p_c2 + p_delta_s2
%
% Because p_delta_s2 is constrained to be nonnegative in the estimation
% script, p_s2 >= p_c2.

theta     = x(1);
theta_dot = x(2);

dx = zeros(2,1);

dx(1) = theta_dot;

% Numerical safety
v_s2  = max(v_s2,  1e-8);
eps_v = max(eps_v, 1e-8);

stribeck_multiplier = exp(-(theta_dot / v_s2)^2);

dry_friction = ...
    (p_c2 + p_delta_s2 * stribeck_multiplier) ...
    * tanh(theta_dot / eps_v);

viscous_friction = p_b2 * theta_dot;

dx(2) = -viscous_friction ...
        -dry_friction ...
        -p_g2*sin(theta);

y = theta;

end