function y = rotpendulumEkfMeasurement(x, ekfInput)
%ROTPENDULUMEKFMEASUREMENT Measurement function for measured-coordinate state.
%
% Explicit Simulink EKF-block signature:
%
%   y = rotpendulumEkfMeasurement(x, ekfInput)
%
% The measurement does not use ekfInput, but the argument is deliberately
% present so the EKF block exposes the same external input vector used by
% the state transition function.
%
% State:
%   x = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]
%
% Output:
%   y = [theta1_meas; theta2_meas]

%#codegen

x = double(x(:));
% ekfInput is intentionally unused. %#ok<NASGU>
y = [x(1); x(2)];

end
