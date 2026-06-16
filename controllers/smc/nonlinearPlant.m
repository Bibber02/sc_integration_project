function xdot = nonlinearPlant(x,u)

% State vector
theta1  = x(1);
theta2  = x(2);
dtheta1 = x(3);
dtheta2 = x(4);

% ============================================================
% Parameters
% ============================================================

pa       = 64.7915;
pb1      = 3.2117e3;
pc1      = 635.3865;
pg1      = 430.5226;
pu       = 2.2424e4;
p0       = 0;

pb2      = 0.0670;
pg2      = 112.0268;
pc2      = 0.2491;
psdelta2 = 0.1836;

epsv1    = 0.05;
epsv2    = 0.03;
vs2      = 3.5;

l1 = 0.1;
g  = 9.81;
pc = (l1/g)*pg2;

% ============================================================
% Dynamics
% ============================================================

c2 = cos(theta2);
s2 = sin(theta2);

% Inertia matrix
M11 = pa + 1 + 2*pc*c2;
M12 = 1 + pc*c2;
M21 = M12;
M22 = 1;

% Coriolis/Centrifugal terms
Cqdot1 = -pc*s2*(2*dtheta1*dtheta2 + dtheta2^2);
Cqdot2 =  pc*s2*dtheta1^2;

% Friction
F1 = pb1*dtheta1 + pc1*tanh(dtheta1/epsv1);

stribeck2 = pc2 + psdelta2*exp(-(dtheta2/vs2)^2);
F2 = pb2*dtheta2 + stribeck2*tanh(dtheta2/epsv2);

% Gravity
G1 = -pg1*sin(theta1) - pg2*sin(theta1 + theta2);
G2 = -pg2*sin(theta1 + theta2);

% Input torques
I1 = pu*u + p0;
I2 = 0;

% Right-hand side
rhs1 = I1 - Cqdot1 - F1 - G1;
rhs2 = I2 - Cqdot2 - F2 - G2;

% Solve M*qdd = rhs
detM = M11*M22 - M12*M21;

ddtheta1 = ( M22*rhs1 - M12*rhs2)/detM;
ddtheta2 = (-M21*rhs1 + M11*rhs2)/detM;

% State derivative
xdot = zeros(4,1);

xdot(1) = dtheta1;
xdot(2) = dtheta2;
xdot(3) = ddtheta1;
xdot(4) = ddtheta2;

end