function [dx, y] = passive_rod_coulomb_model(t, x, u, p_b2, p_g2, p_c2, eps_v, varargin)
%PASSIVE_ROD_COULOMB_MODEL Passive rod model with Coulomb-viscous friction.
% theta_ddot = -p_b2 theta_dot - p_g2 sin(theta) - p_c2 tanh(theta_dot/eps_v)

psi = x(1);
theta_dot = x(2);

eps_v = max(eps_v, 1e-5);

dx = zeros(2, 1);
dx(1) = theta_dot;
dx(2) = -p_b2 * theta_dot ...
        -p_g2 * sin(psi) ...
        -p_c2 * tanh(theta_dot / eps_v);

y = psi;

end
