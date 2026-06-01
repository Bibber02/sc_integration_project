C_filter = [1 0 0 0;
            0 0 1 0];      % measures only theta1 and theta2
D_filter = zeros(2,1);

sys_kf = ss(linsys1.A, linsys1.B, C_filter, D_filter);