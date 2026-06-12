function x_k_1 = rotpendulumEkfStateTransitionMuxed(x_k, ekfInput)

    ekfInput = ekfInput(:);

    u_k = ekfInput(1);
    Ts  = ekfInput(2);
    p   = ekfInput(3:end);

    x_k_1 = rotpendulumEkfStateTransition(x_k, u_k, Ts, p);

end