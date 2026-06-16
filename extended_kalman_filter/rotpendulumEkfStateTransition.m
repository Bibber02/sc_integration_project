function x_k_1 = rotpendulumEkfStateTransition(x_k, varargin)
%ROTPENDULUMEKFSTATETRANSITION Discrete EKF transition using RK4.
%
% Supported call forms:
%
%   x_k_1 = rotpendulumEkfStateTransition(x_k, ekfInput)
%       ekfInput = [u_k; Ts; p(:)]
%
%   x_k_1 = rotpendulumEkfStateTransition(x_k, u_k, Ts, p)
%
% Parameter vector p uses the selected no-p0 down-line convention:
%   [p_a; p_b1; p_c1; p_g1; p_u; ...
%    p_b2; p_g2; p_c2; p_sdelta2; v_s2; eps_v1; eps_v2; ...
%    theta_scale; theta1_offset; theta_abs_down_raw]

x_k = double(x_k(:));

if nargin == 2
    ekfInput = double(varargin{1}(:));
    u_k = ekfInput(1);
    Ts  = ekfInput(2);
    p   = ekfInput(3:end);
elseif nargin == 4
    u_k = double(varargin{1});
    Ts  = double(varargin{2});
    p   = double(varargin{3}(:));
else
    error('rotpendulumEkfStateTransition expects either (x, ekfInput) or (x, u, Ts, p).');
end

[k1, ~] = nonlinearPlant(x_k, u_k, p);
k1 = double(k1(:));

[k2, ~] = nonlinearPlant(x_k + 0.5 * Ts * k1, u_k, p);
k2 = double(k2(:));

[k3, ~] = nonlinearPlant(x_k + 0.5 * Ts * k2, u_k, p);
k3 = double(k3(:));

[k4, ~] = nonlinearPlant(x_k + Ts * k3, u_k, p);
k4 = double(k4(:));

x_k_1 = x_k + (Ts / 6) * (k1 + 2*k2 + 2*k3 + k4);

% EKF block wants exactly 4x1 double.
x_k_1 = double(x_k_1(:));

end
