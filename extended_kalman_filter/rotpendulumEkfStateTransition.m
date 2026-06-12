function x_k_1 = rotpendulumEkfStateTransition(x_k, u_k, Ts, p)

    x_k = double(x_k(:));
    u_k = double(u_k);
    Ts  = double(Ts);
    p   = double(p(:));

    [k1, ~] = nonlinearPlant(x_k, u_k, p);
    k1 = double(k1(:));

    [k2, ~] = nonlinearPlant(x_k + 0.5 * Ts * k1, u_k, p);
    k2 = double(k2(:));

    [k3, ~] = nonlinearPlant(x_k + 0.5 * Ts * k2, u_k, p);
    k3 = double(k3(:));

    [k4, ~] = nonlinearPlant(x_k + Ts * k3, u_k, p);
    k4 = double(k4(:));

    x_k_1 = x_k + (Ts / 6) * (k1 + 2*k2 + 2*k3 + k4);
    x_k_1 = double(x_k_1(:));

end