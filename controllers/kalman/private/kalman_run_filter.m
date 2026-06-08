function [xhatLog, innovationLog] = kalman_run_filter(sysd, Q, R, P0, xhat0, u, y)
%KALMAN_RUN_FILTER Run a covariance-recursive discrete Kalman filter.

A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

n = size(y, 1);
nx = size(A, 1);
ny = size(C, 1);

xhat = xhat0(:);
P = P0;
I = eye(nx);

xhatLog = zeros(n, nx);
innovationLog = zeros(n, ny);

for k = 1:n
    uk = u(k, :).';
    yk = y(k, :).';

    innovation = yk - (C * xhat + D * uk);
    S = C * P * C.' + R;
    K = (P * C.') / S;

    xhat = xhat + K * innovation;
    P = (I - K * C) * P * (I - K * C).' + K * R * K.';
    P = (P + P.') / 2;

    xhatLog(k, :) = xhat.';
    innovationLog(k, :) = innovation.';

    if k < n
        xhat = A * xhat + B * uk;
        P = A * P * A.' + Q;
        P = (P + P.') / 2;
    end
end
end

