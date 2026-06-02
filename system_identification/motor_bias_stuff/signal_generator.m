
Ts = 0.01;                 % Same as your sensor/sample time
A_values = [0.25 0.35 0.45 0.55];
f_tri = 0.25;              % Hz, one period every 20 s

zeroTime = 2;             % seconds between blocks
blockTime = 10;            % seconds per amplitude block

t_all = [];
u_all = [];

t_current = 0;

% Initial zero section
t_block = (0:Ts:zeroTime-Ts).';
u_block = zeros(size(t_block));

t_all = [t_all; t_current + t_block];
u_all = [u_all; u_block];
t_current = t_all(end) + Ts;

for k = 1:length(A_values)

    A = A_values(k);

    % Triangle wave section
    t_block = (0:Ts:blockTime-Ts).';

    % Triangle between -1 and +1
    tri = sawtooth(2*pi*f_tri*t_block, 0.5);

    u_block = A * tri;

    t_all = [t_all; t_current + t_block];
    u_all = [u_all; u_block];
    t_current = t_all(end) + Ts;

    % Zero section after each amplitude
    t_block = (0:Ts:zeroTime-Ts).';
    u_block = zeros(size(t_block));

    t_all = [t_all; t_current + t_block];
    u_all = [u_all; u_block];
    t_current = t_all(end) + Ts;
end

u_motor_bias = timeseries(u_all, t_all);

figure;
plot(t_all, u_all);
grid on;
xlabel('Time [s]');
ylabel('Input command');
title('Motor bias / deadzone / hysteresis identification input');