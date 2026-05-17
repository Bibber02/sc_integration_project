function [dx, y] = pendulum_rod_2_model_lowvel_damping( ...
    t, x, u, p_b2, p_g2, p_c2, p_b_low, v_b, eps_v, varargin)
%PENDULUM_ROD_2_MODEL_LOWVEL_DAMPING
% Passive-link model with Coulomb friction plus low-velocity viscous damping.
%
% Model:
%
%   theta_ddot =
%       -[p_b2 + p_b_low*exp(-(theta_dot/v_b)^2)]*theta_dot
%       -p_c2*tanh(theta_dot/eps_v)
%       -p_g2*sin(theta)
%
% States:
%   x(1) = theta      [rad]
%   x(2) = theta_dot  [rad/s]
%
% Parameters:
%   p_b2     = baseline viscous damping             [1/s]
%   p_g2     = gravity parameter                    [1/s^2]
%   p_c2     = Coulomb friction strength            [rad/s^2]
%   p_b_low  = additional low-velocity damping      [1/s]
%   v_b      = low-velocity damping scale           [rad/s]
%   eps_v    = Coulomb smoothing velocity           [rad/s]
%
% Input:
%   u is included for idnlgrey compatibility, but is not used.

theta     = x(1);
theta_dot = x(2);

dx = zeros(2,1);

dx(1) = theta_dot;

v_b   = max(v_b,   1e-8);
eps_v = max(eps_v, 1e-8);

low_velocity_multiplier = exp(-(theta_dot / v_b)^2);

effective_viscous_damping = ...
    p_b2 + p_b_low * low_velocity_multiplier;

viscous_friction = effective_viscous_damping * theta_dot;
coulomb_friction = p_c2 * tanh(theta_dot / eps_v);

dx(2) = -viscous_friction ...
        -coulomb_friction ...
        -p_g2*sin(theta);

y = theta;

end