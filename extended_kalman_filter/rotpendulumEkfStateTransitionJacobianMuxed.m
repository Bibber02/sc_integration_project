function F = rotpendulumEkfStateTransitionJacobianMuxed(x, ekfInput)

% Discrete-time EKF state-transition Jacobian.
%
% State:
%   x = [theta1; theta2; dtheta1; dtheta2]
%
% Extra EKF input:
%   ekfInput = [u; Ts; p]
%
% The transition function is:
%   x_next = x + Ts * nonlinearPlant(x, u, p)
%
% Therefore:
%   F = d(x_next)/dx = I + Ts * d(f)/dx

u = ekfInput(1);
Ts = ekfInput(2);
p = ekfInput(3:end);

theta1  = x(1);
theta2  = x(2);
dtheta1 = x(3);
dtheta2 = x(4);

pa       = p(1);
pb1      = p(2);
pc1      = p(3);
pg1      = p(4);
% pu     = p(5);
% p0     = p(6);
pb2      = p(7);
pg2      = p(8);
pc2      = p(9);
psdelta2 = p(10);
vs2      = p(11);
epsv1    = p(12);
epsv2    = p(13);
pc       = p(14);

c2 = cos(theta2);
s2 = sin(theta2);

s12 = sin(theta1 + theta2);
c12 = cos(theta1 + theta2);

% Mass matrix terms
M11 = pa + 1 + 2*pc*c2;
M12 = 1 + pc*c2;
M21 = M12;
M22 = 1;

detM = M11*M22 - M12*M21;

% Coriolis/centrifugal terms
velTerm1 = 2*dtheta1*dtheta2 + dtheta2^2;

Cqdot1 = -pc*s2*velTerm1;
Cqdot2 =  pc*s2*dtheta1^2;

% Friction terms
tanh1 = tanh(dtheta1/epsv1);
sech1sq = 1 - tanh1^2;

F1 = pb1*dtheta1 + pc1*tanh1;

exp2 = exp(-(dtheta2/vs2)^2);
stribeck2 = pc2 + psdelta2*exp2;

tanh2 = tanh(dtheta2/epsv2);
sech2sq = 1 - tanh2^2;

F2 = pb2*dtheta2 + stribeck2*tanh2;

% Gravity terms
G1 = -pg1*sin(theta1) - pg2*sin(theta1 + theta2);
G2 = -pg2*sin(theta1 + theta2);

% Input terms
% I1 = pu*u + p0;
% I2 = 0;
% The input does not affect d(f)/dx directly, so it is not needed below.
I1 = p(5)*ekfInput(1) + p(6);
I2 = 0;

rhs1 = I1 - Cqdot1 - F1 - G1;
rhs2 = I2 - Cqdot2 - F2 - G2;

num1 =  M22*rhs1 - M12*rhs2;
num2 = -M21*rhs1 + M11*rhs2;

ddtheta1 = num1 / detM;
ddtheta2 = num2 / detM;

% Derivatives of mass matrix terms
% with respect to [theta1, theta2, dtheta1, dtheta2]
dM11 = zeros(1,4);
dM12 = zeros(1,4);

dM11(2) = -2*pc*s2;
dM12(2) = -pc*s2;

dDetM = zeros(1,4);
dDetM(2) = dM11(2)*M22 - dM12(2)*M21 - M12*dM12(2);

% Since M21 = M12 and M22 = 1:
% dDetM(2) = dM11(2) - 2*M12*dM12(2)

% Derivatives of rhs1
% rhs1 = I1 - Cqdot1 - F1 - G1
drhs1 = zeros(1,4);

% d/dtheta1
drhs1(1) = pg1*cos(theta1) + pg2*c12;

% d/dtheta2
dCqdot1_dtheta2 = -pc*c2*velTerm1;
dG1_dtheta2 = -pg2*c12;
drhs1(2) = -dCqdot1_dtheta2 - dG1_dtheta2;

% d/ddtheta1
dCqdot1_ddtheta1 = -pc*s2*(2*dtheta2);
dF1_ddtheta1 = pb1 + pc1*(1/epsv1)*sech1sq;
drhs1(3) = -dCqdot1_ddtheta1 - dF1_ddtheta1;

% d/ddtheta2
dCqdot1_ddtheta2 = -pc*s2*(2*dtheta1 + 2*dtheta2);
drhs1(4) = -dCqdot1_ddtheta2;

% Derivatives of rhs2
% rhs2 = I2 - Cqdot2 - F2 - G2
drhs2 = zeros(1,4);

% d/dtheta1
drhs2(1) = pg2*c12;

% d/dtheta2
dCqdot2_dtheta2 = pc*c2*dtheta1^2;
dG2_dtheta2 = -pg2*c12;
drhs2(2) = -dCqdot2_dtheta2 - dG2_dtheta2;

% d/ddtheta1
dCqdot2_ddtheta1 = 2*pc*s2*dtheta1;
drhs2(3) = -dCqdot2_ddtheta1;

% d/ddtheta2
dstribeck2_ddtheta2 = psdelta2*exp2*(-2*dtheta2/(vs2^2));
dtanh2_ddtheta2 = (1/epsv2)*sech2sq;
dF2_ddtheta2 = pb2 + dstribeck2_ddtheta2*tanh2 + ...
    stribeck2*dtanh2_ddtheta2;

drhs2(4) = -dF2_ddtheta2;

% Derivatives of acceleration terms
%
% ddtheta1 = num1 / detM
% num1 = rhs1 - M12*rhs2
%
% ddtheta2 = num2 / detM
% num2 = -M12*rhs1 + M11*rhs2
Acont = zeros(4,4);

Acont(1,3) = 1;
Acont(2,4) = 1;

for j = 1:4

    dnum1 = drhs1(j) - dM12(j)*rhs2 - M12*drhs2(j);

    dnum2 = -dM12(j)*rhs1 - M12*drhs1(j) + ...
             dM11(j)*rhs2 + M11*drhs2(j);

    Acont(3,j) = (dnum1*detM - num1*dDetM(j)) / detM^2;
    Acont(4,j) = (dnum2*detM - num2*dDetM(j)) / detM^2;

end

% Discrete-time transition Jacobian
F = eye(4) + Ts*Acont;

end