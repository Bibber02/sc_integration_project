function [dx, y] = passive_rod_viscous_model(t, x, u, p_b2, p_g2, varargin)
%PASSIVE_ROD_VISCOUS_MODEL Passive rod model with viscous damping.
% theta_ddot = -p_b2 theta_dot - p_g2 sin(theta)

psi = x(1);
theta_dot = x(2);

dx = zeros(2, 1);
dx(1) = theta_dot;
dx(2) = -p_b2 * theta_dot - p_g2 * sin(psi);

y = psi;

end
