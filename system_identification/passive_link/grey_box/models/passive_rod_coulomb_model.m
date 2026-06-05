function [dx, y] = passive_rod_coulomb_model(t, x, u, p_b2, p_g2, p_c2, eps_v2, varargin)
%PASSIVE_ROD_COULOMB_MODEL Passive-link free-decay model with Coulomb-viscous friction.
%
% Model:
%   theta2_ddot = -p_b2*theta2_dot
%                 -p_c2*tanh(theta2_dot/eps_v2)
%                 -p_g2*sin(theta2)

%#ok<*INUSD>
theta2 = x(1);
theta2_dot = x(2);

eps_v2 = max(eps_v2, 1e-5);

dx = zeros(2, 1);
dx(1) = theta2_dot;
dx(2) = -p_b2 * theta2_dot ...
        -p_c2 * tanh(theta2_dot / eps_v2) ...
        -p_g2 * sin(theta2);

y = theta2;

end
