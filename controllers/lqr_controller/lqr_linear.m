A = linsys1.A;
B = linsys1.B;

x0 = [pi; 0; 0; 0];
u0 = 0;

% Scaling values: choose expected maximum deviations
theta1_max    = deg2rad(10);
theta1dot_max = 1;
theta2_max    = deg2rad(10);
theta2dot_max = 1;

S = diag([theta1_max, theta1dot_max, theta2_max, theta2dot_max]);

% Scaled linear system
Az = S\A*S;
Bz = S\B;

% Check scaled controllability
Coz = ctrb(Az,Bz);

format long e

disp('Scaled controllability singular values:')
svd(Coz)

disp('Scaled controllability rank:')
rank(Coz)

disp('Scaled controllability condition number:')
cond(Coz)

% LQR
Qz = diag([1, 2, 1, 5]);
R  = 0.1;

Kz = lqr(Az, Bz, Qz, R);

% Convert back to physical-state feedback gain
K = Kz / S;

disp('Physical LQR gain K:')
K

% Closed-loop system in original coordinates
Acl = A - B*K;

disp('Closed-loop poles:')
eig(Acl)