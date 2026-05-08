u = u_Chrip_Amp0_5;
theta1 = theta_1_Chrip_Amp0_5;
theta2 = theta_2_Chrip_Amp0_5;


% figure;
% subplot(3,1,1); plot(u);      ylabel('u');
% subplot(3,1,2); plot(theta1');      ylabel('\theta 1');
% subplot(3,1,3); plot(theta2');      ylabel('\theta 2');
% xlabel('Time [s]')

N = length(u.Time);
h = max(u.Time) / N;
fs = 1/h;

U = fft(u.Data) / N;
Y = fft(theta2.Data) / N;

f = (0:N-1)' * (fs/N);
idx = 1:floor(N/2);

f_plot = f(idx);
G_hat = Y(idx) ./ U(idx);
f_plot

figure;
subplot(2,1,1);
semilogx(20*log10(abs(U)));
ylabel('Magnitude [dB]');
grid on;

subplot(2, 1, 2);
semilogx(angle(G_hat) * 180/pi);
ylabel('Phase [deg]')
xlabel('Frequency [Hz]')
grid on;