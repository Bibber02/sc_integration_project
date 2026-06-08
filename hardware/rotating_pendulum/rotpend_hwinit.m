function cfg = rotpend_hwinit(sampleTime)
%ROTPEND_HWINIT Shared DSCS FPGA I/O conversion constants.

if nargin < 1 || isempty(sampleTime)
    sampleTime = 0.01;
end

cfg.h = sampleTime;

cfg.daoutoffs = [0.00];
cfg.daoutgain = 1 * [-6];

a_theta_1 = 1.200374077783257;
b_theta_1 = 1.070535559513546;

a_theta_2 = 1.222725978474668;
b_theta_2 = 1.190099998448171;

adinoffs = -[b_theta_1 b_theta_2];
adingain = [a_theta_1 a_theta_2];

cfg.adinoffs = [adinoffs 0 0 0 0 0];
cfg.adingain = [adingain 1 1 1 1 1];
end
