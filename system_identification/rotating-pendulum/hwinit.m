%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSCS FPGA interface board: init and I/O conversions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 0.01;
% gains and offsets
daoutoffs = [0.00];                   % output offset
daoutgain = 1*[-6];                   % output gain

a_theta_1 = 1.200374077783257;
b_theta_1 = 1.070535559513546;

a_theta_2 = 1.222725978474668;
b_theta_2 = 1.190099998448171;

% Sensor calibration:
adinoffs = -[b_theta_1 b_theta_2];
adingain = [a_theta_1 a_theta_2];

adinoffs = [adinoffs 0 0 0 0 0];    % input offset
adingain = [adingain 1 1 1 1 1];     % input gain (to radians)
