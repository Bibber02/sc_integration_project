function y = rotpendulumEkfMeasurement(x)
%ROTPENDULUMEKFMEASUREMENT Measurement function for measured-coordinate state.
%
% Simulink EKF-block signature:
%   y = rotpendulumEkfMeasurement(x)
%
% This measurement model does not depend on the control input, Ts, or model
% parameters. The measured outputs are the first two EKF states.
%
% State:
%   x = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]
%
% Output:
%   y = [theta1_meas; theta2_meas]

%#codegen

x = double(x(:));

y = zeros(2,1);
y(1) = x(1);
y(2) = x(2);

y = double(y(:));

end
