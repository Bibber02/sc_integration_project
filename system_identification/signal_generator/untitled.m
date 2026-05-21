Ts        = 0.01;     % sample time [s]
Tsim      = 120;       % total duration [s]
N         = Tsim/Ts + 1;

amplitude = 0.15;     % signal switches between -0.15 and +0.15
clockSamples = 10;    % hold each PRBS value for 25 samples = 0.25 s

Range = [-amplitude amplitude];
Band  = [0 1/clockSamples];

u = idinput(N, 'prbs', Band, Range);

t = (0:N-1)' * Ts;

u_ts = timeseries(u, t);
u_ts.Name = 'u_PRBS';

figure;
stairs(t, u);
grid on;
xlabel('Time [s]');
ylabel('u');
title('PRBS generated with idinput');