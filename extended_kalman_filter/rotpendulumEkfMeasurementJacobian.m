function H = rotpendulumEkfMeasurementJacobian(x)

% Measurement Jacobian for:
%
%   y = [theta1; theta2]
%
% State:
%   x = [theta1; theta2; dtheta1; dtheta2]

x = x(:);

H = [1 0 0 0;
     0 1 0 0];

end