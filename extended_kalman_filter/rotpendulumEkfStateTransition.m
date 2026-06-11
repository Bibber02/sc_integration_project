function x_k_1 = rotpendulumEkfStateTransition(x_k, u_k, Ts, p)

    x_k = x_k(:);

    [k1, ~] = nonlinearPlant(x_k, u_k, p);
    [k2, ~] = nonlinearPlant(x_k + 0.5 * Ts * k1, u_k, p);
    [k3, ~] = nonlinearPlant(x_k + 0.5 * Ts * k2, u_k, p);
    [k4, ~] = nonlinearPlant(x_k + Ts * k3, u_k, p);

    x_k_1 = x_k + (Ts / 6) * (k1 + 2*k2 + 2*k3 + k4);

end