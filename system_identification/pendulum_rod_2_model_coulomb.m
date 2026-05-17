function [dx, y] = pendulum_rod_2_model_coulomb( ...
    t, x, u, p_b2, p_g2, p_c2, eps_v, varargin)
%PENDULUM_ROD_2_MODEL_COULOMB Nonlinear passive-link model for idnlgrey.
%
% Model:
%
%   theta_ddot =
%       -p_b2*theta_dot
%       -p_c2*tanh(theta_dot/eps_v)
%       -p_g2*sin(theta)
%
% States:
%   x(1) = theta      [rad]
%   x(2) = theta_dot  [rad/s]
%
% Parameters:
%   p_b2  = viscous damping parameter         [1/s]
%   p_g2  = gravity parameter                 [1/s^2]
%   p_c2  = smooth Coulomb friction strength  [rad/s^2]
%   eps_v = smoothing velocity                [rad/s]
%
% Input:
%   u is included for idnlgrey compatibility, but is not used.

theta     = x(1);
theta_dot = x(2);

dx = zeros(2,1);

dx(1) = theta_dot;

% Numerical safety
eps_v = max(eps_v, 1e-8);

viscous_friction = p_b2 * theta_dot;
coulomb_friction = p_c2 * tanh(theta_dot / eps_v);

dx(2) = -viscous_friction ...
        -coulomb_friction ...
        -p_g2*sin(theta);

y = theta;

end