h = 0.01;
Tsim = 90;
N = Tsim / h;

u_prbs = idinput(N, 'prbs', [3*h, 20*h], [-1, 1]);


t = (0:h:Tsim-h)';
simin = [t, u_prbs];