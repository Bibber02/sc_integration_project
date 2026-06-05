% Extract discrete matrices from the LTI system
Ad = sys_disc.A;
Bd = sys_disc.B;
Cd = sys_disc.C;
Dd = sys_disc.D;
Ts = sys_disc.Ts;

% State order:
% dx = [dtheta1; dtheta2; dtheta1_dot; dtheta2_dot]

Q_lqr = diag([
    5;    % theta1 error
    1;    % theta2 error
    0.1;      % theta1_dot
    0.1       % theta2_dot
]);

R_lqr = 5;

% Discrete-time LQR gain
K_lqr = dlqr(Ad, Bd, Q_lqr, R_lqr);

disp('K_lqr = ');
disp(K_lqr);

disp('Closed-loop poles:');
disp(eig(Ad - Bd*K_lqr));

fs = 1/h;
fc = 10;
wn = fc / (fs/2);
[b, a] = butter(1, wn);