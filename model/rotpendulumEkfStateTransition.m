function xNext = rotpendulumEkfStateTransition(x, ekfInput)
%#codegen
%ROTPENDULUMEKFSTATETRANSITION One-step RK4 transition for the EKF.

uDev = ekfInput(1);
Ts = ekfInput(2);
u0 = ekfInput(3);
xOffset = ekfInput(4:7);
p = ekfInput(8:21);

xAbs = x(:) + xOffset(:);
uAbs = u0 + uDev;

k1 = nonlinearPlant(xAbs, uAbs, p);
k2 = nonlinearPlant(xAbs + 0.5 * Ts * k1, uAbs, p);
k3 = nonlinearPlant(xAbs + 0.5 * Ts * k2, uAbs, p);
k4 = nonlinearPlant(xAbs + Ts * k3, uAbs, p);

xNextAbs = xAbs + (Ts / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
xNext = xNextAbs - xOffset(:);
end
