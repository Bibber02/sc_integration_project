function x_k_1 = rotpendulumEkfStateTransitionMuxed(x_k, ekfInput)
%ROTPENDULUMEKFSTATETRANSITIONMUXED Wrapper for Simulink EKF block.
%
% ekfInput = [u; Ts; p(:)]

x_k = double(x_k(:));
ekfInput = double(ekfInput(:));

x_k_1 = rotpendulumEkfStateTransition(x_k, ekfInput);
x_k_1 = double(x_k_1(:));

end
