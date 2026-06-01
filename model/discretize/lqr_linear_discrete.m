linsys1_d = c2d(linsys1, h, 'zoh');  % Zero-Order Hold

% Extract discrete matrices
Ad = linsys1_d.A;
Bd = linsys1_d.B;
Cd = linsys1_d.C;

Qz_d = diag([3, 1, 5, 1]);
R_d  = 1;

% Redesign LQR in discrete time
Kx = dlqr(Ad, Bd, Qz_d, R_d);


C_filter = [1 0 0 0;
            0 0 1 0];  
D_filter = zeros(2,1);

sys_kf_c = ss(linsys1.A, linsys1.B, C_filter, D_filter);

% Discretize
sys_kf_d = c2d(sys_kf_c, h, 'zoh');

K
Kx