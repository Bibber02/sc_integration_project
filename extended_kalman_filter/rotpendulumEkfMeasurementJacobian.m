function H = rotpendulumEkfMeasurementJacobian(x, ekfInput)
%ROTPENDULUMEKFMEASUREMENTJACOBIAN Measurement Jacobian.
%
% Explicit Simulink EKF-block signature:
%   H = rotpendulumEkfMeasurementJacobian(x, ekfInput)
%
% ekfInput is intentionally unused, but included to match the EKF block's
% external input argument convention.
%
% Measurement:
%   y = [theta1_meas; theta2_meas]
%
% State:
%   x = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]

%#codegen

x = x(:); %#ok<NASGU>
% ekfInput intentionally unused. %#ok<NASGU>

H = [1 0 0 0;
     0 1 0 0];

end
