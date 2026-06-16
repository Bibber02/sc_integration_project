function [xdot, y] = nonlinearPlant(x, u, p)
%NONLINEARPLANT Reduced full-system model using the selected no-p0
%down-line calibrated identification result.
%
% State convention:
%   x = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]
%
% Parameter vector:
%   p = [p_a; p_b1; p_c1; p_g1; p_u; ...
%        p_b2; p_g2; p_c2; p_sdelta2; v_s2; eps_v1; eps_v2; ...
%        theta_scale; theta1_offset; theta_abs_down_raw]
%
% The dynamics are evaluated in corrected physical coordinates, but the
% state and output remain in the measured/saved coordinate. There is no
% motor-bias p0 term.

%#codegen

theta1_m  = x(1);
theta2_m  = x(2);
dtheta1_m = x(3);
dtheta2_m = x(4);

pa       = p(1);
pb1      = p(2);
pc1      = p(3);
pg1      = p(4);
pu       = p(5);
pb2      = p(6);
pg2      = p(7);
pc2      = p(8);
psdelta2 = p(9);
vs2      = p(10);
epsv1    = p(11);
epsv2    = p(12);

theta_scale        = p(13);
theta1_offset      = p(14);
theta_abs_down_raw = p(15);

% Known constants used in the reduced model
g  = 9.81;
l1 = 0.10;

theta_scale = max(theta_scale, 1e-6);
epsv1 = max(epsv1, 1e-5);
epsv2 = max(epsv2, 1e-5);
vs2   = max(vs2,   1e-5);

% Coupling parameter and dependent theta2 offset from the down-line
% calibration constraint.
pc = l1 * pg2 / g;
theta2_offset = pi - theta_scale * theta_abs_down_raw - theta1_offset;

% Corrected physical coordinates
theta1  = theta_scale * theta1_m  + theta1_offset;
theta2  = theta_scale * theta2_m  + theta2_offset;
dtheta1 = theta_scale * dtheta1_m;
dtheta2 = theta_scale * dtheta2_m;

c2 = cos(theta2);
s2 = sin(theta2);

% Reduced inertia matrix
M11 = pa + 1 + 2*pc*c2;
M12 = 1 + pc*c2;
M21 = M12;
M22 = 1;

% Reduced Coriolis/centrifugal vector C(q,qdot)qdot
Cqdot1 = -pc*s2*(2*dtheta1*dtheta2 + dtheta2^2);
Cqdot2 =  pc*s2*dtheta1^2;

% Friction
F1 = pb1*dtheta1 + pc1*tanh(dtheta1/epsv1);

stribeck2 = pc2 + psdelta2*exp(-(dtheta2/vs2)^2);
F2 = pb2*dtheta2 + stribeck2*tanh(dtheta2/epsv2);

% Gravity
G1 = -pg1*sin(theta1) - pg2*sin(theta1 + theta2);
G2 = -pg2*sin(theta1 + theta2);

% Input: no constant motor-bias term
I1 = pu * u;
I2 = 0;

rhs1 = I1 - Cqdot1 - F1 - G1;
rhs2 = I2 - Cqdot2 - F2 - G2;

detM = M11*M22 - M12*M21;

ddtheta1_phys = ( M22*rhs1 - M12*rhs2) / detM;
ddtheta2_phys = (-M21*rhs1 + M11*rhs2) / detM;

% Convert physical accelerations back to measured-coordinate accelerations.
xdot = zeros(4,1);
xdot(1) = dtheta1_m;
xdot(2) = dtheta2_m;
xdot(3) = ddtheta1_phys / theta_scale;
xdot(4) = ddtheta2_phys / theta_scale;

% Sensors/output are in measured/saved coordinates.
y = zeros(2,1);
y(1) = theta1_m;
y(2) = theta2_m;

end
