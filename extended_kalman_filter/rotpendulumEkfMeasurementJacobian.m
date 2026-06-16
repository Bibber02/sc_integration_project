function H = rotpendulumEkfMeasurementJacobian(x)
%ROTPENDULUMEKFMEASUREMENTJACOBIAN Measurement Jacobian.
%
% Simulink EKF-block signature:
%   H = rotpendulumEkfMeasurementJacobian(x)
%
% Measurement:
%   y = [theta1_meas; theta2_meas]
%
% State:
%   x = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]
%
% The measurement does not depend on the control input, Ts, or parameters.

%#codegen

x = x(:); %#ok<NASGU>

H = [1 0 0 0;
     0 1 0 0];

end
