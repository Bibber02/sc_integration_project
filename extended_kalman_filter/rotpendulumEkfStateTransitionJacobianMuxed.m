function F = rotpendulumEkfStateTransitionJacobianMuxed(x, ekfInput)
%ROTPENDULUMEKFSTATETRANSITIONJACOBIANMUXED Discrete-time EKF Jacobian.
%
% This numerical Jacobian is intentionally tied to the actual RK4 transition
% function rotpendulumEkfStateTransition. That makes it consistent with the
% selected no-p0 down-line calibrated nonlinearPlant, including theta_scale
% and theta offsets.
%
% State:
%   x = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]
%
% Extra EKF input:
%   ekfInput = [u; Ts; p(:)]

x = double(x(:));
ekfInput = double(ekfInput(:));

nx = 4;
F = zeros(nx, nx);

baseStep = 1e-6;

for k = 1:nx
    dx = zeros(nx, 1);
    h = baseStep * max(1, abs(x(k)));
    dx(k) = h;

    xPlus  = rotpendulumEkfStateTransition(x + dx, ekfInput);
    xMinus = rotpendulumEkfStateTransition(x - dx, ekfInput);

    F(:, k) = (xPlus - xMinus) / (2*h);
end

end
