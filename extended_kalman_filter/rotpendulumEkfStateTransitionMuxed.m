function x_k_1 = rotpendulumEkfStateTransitionMuxed(x_k, ekfInput)

    x_k = double(x_k(:));
    ekfInput = double(ekfInput(:));

    u_k = ekfInput(1);
    Ts  = ekfInput(2);
    p   = ekfInput(3:end);

    x_k_1 = rotpendulumEkfStateTransition(x_k, u_k, Ts, p);
    x_k_1 = double(x_k_1(:));

end