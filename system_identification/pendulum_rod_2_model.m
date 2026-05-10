function [dx, y] = pendulum_rod_2_model(t, x, u, p_b2, p_g2, varargin)
%PENDULUM_ROD_2_MODEL Nonlinear free pendulum model for idnlgrey.
%
% States:
%   x(1) = theta      [rad]
%   x(2) = theta_dot  [rad/s]
%
% Parameters:
%   p_b2 = damping parameter
%   p_g2 = gravity parameter
%   y0   = measurement offset
%
% Input:
%   u is included for compatibility, but not used.
%
% varargin:
%   included for compatibility with idnlgrey FileArgument.

theta     = x(1);
theta_dot = x(2);

dx = zeros(2,1);

dx(1) = theta_dot;
dx(2) = -p_b2*theta_dot + p_g2*sin(theta);

y = theta;

end