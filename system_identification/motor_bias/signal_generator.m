Ts = 0.01;

% Motor deadband starts around 0.2, so stay above it.
A_values = [0.2 0.25 0.3];

zeroStart = 3;     % initial zero time [s]
holdTime  = 2;     % command hold time [s]
zeroTime  = 2;     % zero between pulses [s]
nRepeats  = 2;     % repeat full amplitude sequence

t_all = [];
u_all = [];
t_current = 0;

% Initial zero
t_block = (0:Ts:zeroStart-Ts).';
u_block = zeros(size(t_block));

t_all = [t_all; t_current + t_block];
u_all = [u_all; u_block];
t_current = t_all(end) + Ts;

for r = 1:nRepeats
    for k = 1:length(A_values)

        A = A_values(k);

        % Positive pulse
        t_block = (0:Ts:holdTime-Ts).';
        u_block = A * ones(size(t_block));

        t_all = [t_all; t_current + t_block];
        u_all = [u_all; u_block];
        t_current = t_all(end) + Ts;

        % Zero
        t_block = (0:Ts:zeroTime-Ts).';
        u_block = zeros(size(t_block));

        t_all = [t_all; t_current + t_block];
        u_all = [u_all; u_block];
        t_current = t_all(end) + Ts;

        % Negative pulse
        t_block = (0:Ts:holdTime-Ts).';
        u_block = -A * ones(size(t_block));

        t_all = [t_all; t_current + t_block];
        u_all = [u_all; u_block];
        t_current = t_all(end) + Ts;

        % Zero
        t_block = (0:Ts:zeroTime-Ts).';
        u_block = zeros(size(t_block));

        t_all = [t_all; t_current + t_block];
        u_all = [u_all; u_block];
        t_current = t_all(end) + Ts;
    end
end

u_motor_bias = timeseries(u_all, t_all);

figure;
plot(t_all, u_all, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Input command');
title('Motor bias ID input: paired positive/negative pulses');