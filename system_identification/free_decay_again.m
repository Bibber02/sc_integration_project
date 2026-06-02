ts_chopped = getsampleusingtime(theta_2, 4, 12);

t = ts_chopped.Time;
y = ts_chopped.Data;

f = fittype('A * exp(-alpha * x) .* sin(omega * x + phi) + b', 'independent', 'x', 'coefficients', {'A', 'alpha', 'omega', 'phi', 'b'});
x0 = [0.3, 0.05, 2*pi*1, 0, 0.01];

fitted = fit(t(:), y(:), f, 'StartPoint', x0);

alpha = fitted.alpha;
tau = 1 / alpha;

disp(fitted)
plot(fitted, t, y)