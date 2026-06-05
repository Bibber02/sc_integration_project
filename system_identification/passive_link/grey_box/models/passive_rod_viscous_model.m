function [dx, y] = passive_rod_viscous_model(t, x, u, p_b2, p_g2, varargin)
%PASSIVE_ROD_VISCOUS_MODEL Passive-link free-decay model with viscous friction.
%
% States:
%   x(1) = theta2       [rad]
%   x(2) = theta2_dot   [rad/s]
%
% Model:
%   theta2_ddot = -p_b2*theta2_dot - p_g2*sin(theta2)

%#ok<*INUSD>
theta2 = x(1);
theta2_dot = x(2);

dx = zeros(2, 1);
dx(1) = theta2_dot;
dx(2) = -p_b2 * theta2_dot - p_g2 * sin(theta2);

y = theta2;

end
