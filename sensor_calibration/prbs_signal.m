h = 0.001;
Tsim = 90;
N = Tsim / h;

u_prbs = idinput(N, 'prbs', [0, 0.05], [-1, 1]);


t = (0:h:Tsim-h)';
simin = [t, u_prbs];