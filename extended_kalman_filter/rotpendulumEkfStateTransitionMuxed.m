function x_k_1 = rotpendulumEkfStateTransitionMuxed(x_k, ekfInput)
%ROTPENDULUMEKFSTATETRANSITIONMUXED Compatibility wrapper.
%
% ekfInput = [u_model; Ts; p(:)]
%
% This wrapper has the same explicit two-input signature as the main EKF
% transition function. It is safe to use in the EKF block, but the preferred
% function name is rotpendulumEkfStateTransition.

%#codegen

x_k_1 = rotpendulumEkfStateTransition(x_k, ekfInput);
x_k_1 = double(x_k_1(:));

end
