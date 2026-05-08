h = 0.005;
Tsim = 90;
N = Tsim / h;
t = (0:h:Tsim-h)';
A = 1;

freqs = logspace(-2, 1.5, 25);
omega = freqs;
u_multi = zeros(length(t), 1);
for k = 1:length(freqs)
    phi = 2*pi*rand();
    u_multi = u_multi + sin(freqs(k)*t+phi);
end
u_multi = A * u_multi / max(abs(u_multi));

simin = [t, u_multi];