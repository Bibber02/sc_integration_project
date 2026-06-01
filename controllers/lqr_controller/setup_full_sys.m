h = 0.01;
hwinit;
build_sysid_parameters;
load('linsys1.mat');
sysd = c2d(linsys1, h, 'zoh')
lqr_linear;
setup_kalman;
sys_kf_d = c2d(sys_kf, h, 'zoh')
open_system('MainPendulum.slx')
