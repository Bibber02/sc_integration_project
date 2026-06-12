function x_k_1 = rotpendulumEkfStateTransition(x_k, ekfInput)
%ROTPENDULUMEKFSTATETRANSITION State transition function for Simulink EKF block.
%
% The Simulink EKF block calls this function as:
%
%   x_k_1 = rotpendulumEkfStateTransition(x_k, ekfInput)
%
% where:
%
%   x_k      = current EKF state estimate, supplied internally by EKF block
%   ekfInput = [u_k; Ts; p(:)]
%
% So the EKF extra input port should receive one vector:
%
%   [u; Ts; p(:)]

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

    % Important for the Simulink EKF block:
    % output must be exactly 4x1 double
    x_k_1 = double(x_k_1(:));

end