function x_k_1 = rotpendulumEkfStateTransition(x_k, u_k, Ts, p)

    [x_dot, ~] = nonlinearPlant(x_k, u_k, p);
    x_k_1 = x_k + x_dot * Ts;

end
