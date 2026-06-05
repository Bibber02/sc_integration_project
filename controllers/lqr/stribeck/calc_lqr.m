% Extract discrete matrices from the LTI system
Ad = sys_disc.A;
Bd = sys_disc.B;
Cd = sys_disc.C;
Dd = sys_disc.D;
Ts = sys_disc.Ts;

% State order:
% dx = [dtheta1; dtheta2; dtheta1_dot; dtheta2_dot]

Q_lqr = diag([
    100;    % theta1 error
    100;    % theta2 error
    1;      % theta1_dot
    1       % theta2_dot
]);

R_lqr = 1;

% Discrete-time LQR gain
K_lqr = dlqr(Ad, Bd, Q_lqr, R_lqr);

disp('K_lqr = ');
disp(K_lqr);

disp('Closed-loop poles:');
disp(eig(Ad - Bd*K_lqr));