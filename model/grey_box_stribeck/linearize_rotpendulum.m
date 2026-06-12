function lin = linearize_rotpendulum_simple(x0, Ts, p, u0)

x0 = x0(:);

[f0, y0] = nonlinearPlant(x0, u0, p);

hx = 1e-6;
hu = 1e-6;

nx = numel(x0);
ny = numel(y0);

A = zeros(nx,nx);
B = zeros(nx,1);
C = zeros(ny,nx);
D = zeros(ny,1);

for i = 1:nx
    dx = zeros(nx,1);
    dx(i) = hx;

    [fp, yp] = nonlinearPlant(x0 + dx, u0, p);
    [fm, ym] = nonlinearPlant(x0 - dx, u0, p);

    A(:,i) = (fp - fm)/(2*hx);
    C(:,i) = (yp - ym)/(2*hx);
end

[fp, yp] = nonlinearPlant(x0, u0 + hu, p);
[fm, ym] = nonlinearPlant(x0, u0 - hu, p);

B(:,1) = (fp - fm)/(2*hu);
D(:,1) = (yp - ym)/(2*hu);

sys_lin = ss(A,B,C,D);
sys_lin.StateName  = {'theta1','theta2','theta1_dot','theta2_dot'};
sys_lin.InputName  = {'u'};
sys_lin.OutputName = {'theta1','theta2'};

sys_disc = c2d(sys_lin, Ts, 'zoh');

lin.A = A;
lin.B = B;
lin.C = C;
lin.D = D;
lin.Ad = sys_disc.A;
lin.Bd = sys_disc.B;
lin.Cd = sys_disc.C;
lin.Dd = sys_disc.D;
lin.sys_lin = sys_lin;
lin.sys_disc = sys_disc;
lin.Ts = Ts;
lin.x0 = x0;
lin.u0 = u0;
lin.f0 = f0;
lin.y0 = y0;

end