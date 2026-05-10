function [dx, y] = pendulum_rod_2_model_coulomb(t, x, u, p_b2, p_g2, p_c2, varargin)
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
%   p_b2 = viscous damping parameter        [1/s]
%   p_g2 = gravity parameter                [1/s^2]
%   p_c2 = smooth Coulomb friction strength [rad/s^2]
%
% Input:
%   u is included for idnlgrey compatibility, but is not used.
%
% Notes:
%   theta is assumed to be the angle deviation from the stable hanging
%   equilibrium. Therefore the gravitational term is -p_g2*sin(theta).
%
%   tanh(theta_dot/eps_v) is used instead of sign(theta_dot) because sign()
%   is discontinuous and can make nonlinear grey-box estimation difficult.

theta     = x(1);
theta_dot = x(2);

% Smooth transition velocity for Coulomb friction.
% Smaller values approximate sign(theta_dot) more closely.
% Larger values make the model smoother and easier to estimate.
eps_v = 0.01;  % [rad/s]

dx = zeros(2,1);

dx(1) = theta_dot;

dx(2) = -p_b2*theta_dot ...
        -p_c2*tanh(theta_dot/eps_v) ...
        -p_g2*sin(theta);

y = theta;

end