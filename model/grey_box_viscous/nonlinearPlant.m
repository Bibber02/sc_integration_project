function [xdot, y] = nonlinearPlant(x, u, p)
%#codegen
theta1  = x(1);
theta2  = x(2);
dtheta1 = x(3);
dtheta2 = x(4);

pa  = p(1);
pb1 = p(2);
pg1 = p(3);
pu  = p(4);
p0  = p(5);
pb2 = p(6);
pg2 = p(7);
pc  = p(8);

c2 = cos(theta2);
s2 = sin(theta2);

M11 = pa + 1 + 2*pc*c2;
M12 = 1 + pc*c2;
M21 = M12;
M22 = 1;

Cqdot1 = -pc*s2*(2*dtheta1*dtheta2 + dtheta2^2);
Cqdot2 =  pc*s2*dtheta1^2;

F1 = pb1*dtheta1;
F2 = pb2*dtheta2;

G1 = -pg1*sin(theta1) - pg2*sin(theta1 + theta2);
G2 = -pg2*sin(theta1 + theta2);

I1 = pu * u + p0;
I2 = 0;

rhs1 = I1 - Cqdot1 - F1 - G1;
rhs2 = I2 - Cqdot2 - F2 - G2;

detM = M11*M22 - M12*M21;

ddtheta1 = ( M22*rhs1 - M12*rhs2) / detM;
ddtheta2 = (-M21*rhs1 + M11*rhs2) / detM;

xdot = zeros(4,1);
xdot(1) = dtheta1;
xdot(2) = dtheta2;
xdot(3) = ddtheta1;
xdot(4) = ddtheta2;

y = zeros(2,1);
y(1) = theta1;
y(2) = theta2;
end
