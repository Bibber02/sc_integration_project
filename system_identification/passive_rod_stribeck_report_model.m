function [dx, y] = passive_rod_stribeck_report_model(t, x, u, p_b2, p_g2, p_c2, p_sdelta, v_s, eps_v, varargin)
%PASSIVE_ROD_STRIBECK_REPORT_MODEL Report-derived passive rod Stribeck model.
% p_sdelta = F_S - F_C, so the static level is p_c2 + p_sdelta.
% F_Str = p_b2*theta_dot + (p_c2 + p_sdelta*exp(-(theta_dot/v_s)^2))*tanh(theta_dot/eps_v)
% theta_ddot = -F_Str - p_g2*sin(theta)

psi = x(1);
theta_dot = x(2);

v_s = max(v_s, 1e-5);
eps_v = max(eps_v, 1e-5);

friction_level = p_c2 + p_sdelta * exp(-(theta_dot / v_s)^2);
friction = p_b2 * theta_dot + friction_level * tanh(theta_dot / eps_v);

dx = zeros(2, 1);
dx(1) = theta_dot;
dx(2) = -friction - p_g2 * sin(psi);

y = psi;

end
