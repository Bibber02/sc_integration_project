function x_k_1 = rotpendulumEkfStateTransition(x_k, ekfInput)
%ROTPENDULUMEKFSTATETRANSITION Discrete EKF transition using RK4.
%
% Explicit Simulink EKF-block signature:
%
%   x_k_1 = rotpendulumEkfStateTransition(x_k, ekfInput)
%
% where:
%
%   x_k      = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]
%   ekfInput = [u_model; Ts; p(:)]
%
% The explicit two-input signature is intentional. Do not use varargin here,
% because the Simulink EKF block infers external input ports from the
% function argument list.
%
% Parameter vector p uses the selected no-p0 down-line convention:
%   [p_a; p_b1; p_c1; p_g1; p_u; ...
%    p_b2; p_g2; p_c2; p_sdelta2; v_s2; eps_v1; eps_v2; ...
%    theta_scale; theta1_offset; theta_abs_down_raw]

%#codegen

x_k = double(x_k(:));
ekfInput = double(ekfInput(:));

u_k = ekfInput(1);
Ts  = ekfInput(2);
p   = ekfInput(3:end);

[k1, ~] = nonlinearPlant(x_k, u_k, p);
k1 = double(k1(:));

[k2, ~] = nonlinearPlant(x_k + 0.5 * Ts * k1, u_k, p);
k2 = double(k2(:));

[k3, ~] = nonlinearPlant(x_k + 0.5 * Ts * k2, u_k, p);
k3 = double(k3(:));

[k4, ~] = nonlinearPlant(x_k + Ts * k3, u_k, p);
k4 = double(k4(:));

x_k_1 = x_k + (Ts / 6) * (k1 + 2*k2 + 2*k3 + k4);
x_k_1 = double(x_k_1(:));

end
