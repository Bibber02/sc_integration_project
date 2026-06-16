function H = rotpendulumEkfMeasurementJacobian(x)
%ROTPENDULUMEKFMEASUREMENTJACOBIAN Measurement Jacobian.
%
% Measurement:
%   y = [theta1_meas; theta2_meas]
%
% State:
%   x = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]

x = x(:); %#ok<NASGU>

H = [1 0 0 0;
     0 1 0 0];

end
