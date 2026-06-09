function y = rotpendulumEkfMeasurement(x)
%#codegen
%ROTPENDULUMEKFMEASUREMENT Angle measurement model for the EKF.

y = [x(1); x(2)];
end
