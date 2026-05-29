Ts        = 0.01;     % sample time [s]
Tsim      = 121;       % total duration [s]
delayTime = 1.0;

N         = Tsim/Ts + 1;


amplitude = 0.15;     % signal switches between -0.15 and +0.15
clockSamples = 10;    % hold each PRBS value for 25 samples = 0.25 s

Range = [-amplitude amplitude];
Band  = [0 1/clockSamples];

N_delay = round(delayTime / Ts);

N_prbs = N - N_delay;

u_prbs = idinput(N_prbs, 'prbs', Band, Range);

u = [zeros(N_delay, 1); u_prbs];

t = (0:N-1)' * Ts;

u_ts = timeseries(u, t);
u_ts.Name = 'u_PRBS';

figure;
stairs(t, u);
grid on;
xlabel('Time [s]');
ylabel('u');
title('PRBS generated with idinput');