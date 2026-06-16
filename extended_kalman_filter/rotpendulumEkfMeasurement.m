function y = rotpendulumEkfMeasurement(x)
%ROTPENDULUMEKFMEASUREMENT Measurement function for measured-coordinate state.
%
% State:
%   x = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]
%
% Output:
%   y = [theta1_meas; theta2_meas]

x = x(:);
y = [x(1); x(2)];
end
