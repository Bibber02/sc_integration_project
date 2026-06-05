clear;
clc;
close all;

% ------------------------------------------------------------
% 1. Load parameter vector p
% ------------------------------------------------------------
% This creates p in the workspace.
%
% Viscous full-system parameter vector:
%   p = [p_a; p_b1; p_g1; p_u; p_0; p_b2; p_g2; p_c]
%
% This is the viscous-only counterpart of the earlier Stribeck linearization
% script. The dynamics are the same reduced full-system dynamics, but the
% friction terms are:
%   F1 = p_b1 * theta1_dot
%   F2 = p_b2 * theta2_dot

load_parameters;

% ------------------------------------------------------------
% 2. Define operating point
% ------------------------------------------------------------
% State order:
% x = [theta1; theta2; theta1_dot; theta2_dot]
%
% This keeps the same default operating point as the Stribeck linearization
% script. Change x0 here if you use a different angle convention.

x0 = [0; 0; 0; 0];

% Input equation:
% I1 = p_u*u + p_0
%
% Therefore, to cancel the torque offset at equilibrium:
% p_u*u0 + p_0 = 0
%
% p(4) = p_u
% p(5) = p_0

u0 = -p(5)/p(4);

% ------------------------------------------------------------
% 3. Check nonlinear model at operating point
% ------------------------------------------------------------

[f0, y0] = nonlinearPlantViscous(x0, u0, p);

disp('Operating point check:');
disp('x0 = ');
disp(x0);

disp('u0 = ');
disp(u0);

disp('f(x0,u0) = ');
disp(f0);

disp('y0 = ');
disp(y0);

if norm(f0) > 1e-6
    warning('The selected operating point is not an exact equilibrium. The linear model will have an affine offset.');
end

% ------------------------------------------------------------
% 4. Numerical linearization settings
% ------------------------------------------------------------

nx = numel(x0);
nu = 1;
ny = numel(y0);

A = zeros(nx,nx);
B = zeros(nx,nu);
C = zeros(ny,nx);
D = zeros(ny,nu);

% Finite-difference step sizes.
hx = 1e-6;
hu = 1e-6;

% ------------------------------------------------------------
% 5. Compute A and C matrices
% ------------------------------------------------------------
% A = df/dx
% C = dy/dx

for i = 1:nx
    dx = zeros(nx,1);
    dx(i) = hx;

    [f_plus,  y_plus]  = nonlinearPlantViscous(x0 + dx, u0, p);
    [f_minus, y_minus] = nonlinearPlantViscous(x0 - dx, u0, p);

    A(:,i) = (f_plus - f_minus)/(2*hx);
    C(:,i) = (y_plus - y_minus)/(2*hx);
end

% ------------------------------------------------------------
% 6. Compute B and D matrices
% ------------------------------------------------------------
% B = df/du
% D = dy/du

[f_plus,  y_plus]  = nonlinearPlantViscous(x0, u0 + hu, p);
[f_minus, y_minus] = nonlinearPlantViscous(x0, u0 - hu, p);

B(:,1) = (f_plus - f_minus)/(2*hu);
D(:,1) = (y_plus - y_minus)/(2*hu);

% ------------------------------------------------------------
% 7. Create continuous-time state-space model
% ------------------------------------------------------------

sys_lin = ss(A,B,C,D);

sys_lin.StateName = {'theta1','theta2','theta1_dot','theta2_dot'};
sys_lin.InputName = {'u'};
sys_lin.OutputName = {'theta1','theta2'};

% ------------------------------------------------------------
% 8. Discretize continuous-time linear model
% ------------------------------------------------------------

Ts = 0.01;   % Sampling period [s]

% Zero-order hold discretization.
% This assumes the input u is held constant during each sample interval.
sys_disc = c2d(sys_lin, Ts, 'zoh');

Ad = sys_disc.A;
Bd = sys_disc.B;
Cd = sys_disc.C;
Dd = sys_disc.D;

% ------------------------------------------------------------
% 9. Display results
% ------------------------------------------------------------

disp('Continuous-time viscous linear model:');
disp('A = ');
disp(A);

disp('B = ');
disp(B);

disp('C = ');
disp(C);

disp('D = ');
disp(D);

disp('Eigenvalues of continuous-time A:');
disp(eig(A));

disp(' ');
disp('Discrete-time viscous linear model:');
disp('Ad = ');
disp(Ad);

disp('Bd = ');
disp(Bd);

disp('Cd = ');
disp(Cd);

disp('Dd = ');
disp(Dd);

disp('Eigenvalues of discrete-time Ad:');
disp(eig(Ad));

% ------------------------------------------------------------
% 10. Save result
% ------------------------------------------------------------

save('linearized_plant_viscous.mat', ...
     'A','B','C','D', ...
     'Ad','Bd','Cd','Dd', ...
     'sys_lin','sys_disc', ...
     'Ts','x0','u0','p');

disp('Saved continuous and discrete viscous linearized models to linearized_plant_viscous.mat');

function [xdot, y] = nonlinearPlantViscous(x, u, p)
% nonlinearPlantViscous
%
% Reduced nonlinear plant model for the rotational double pendulum using
% viscous friction only.
%
% Inputs:
%   x = [theta1; theta2; theta1_dot; theta2_dot]
%   u = motor input command / voltage-like input
%   p = [p_a; p_b1; p_g1; p_u; p_0; p_b2; p_g2; p_c]
%
% Outputs:
%   xdot = state derivative
%   y    = measured output [theta1; theta2]

% ------------------------------------------------------------
% 1. States
% ------------------------------------------------------------
theta1  = x(1);
theta2  = x(2);
dtheta1 = x(3);
dtheta2 = x(4);

% ------------------------------------------------------------
% 2. Unpack parameter vector
% ------------------------------------------------------------
p_a  = p(1);

p_b1 = p(2);
p_g1 = p(3);
p_u  = p(4);
p_0  = p(5);

p_b2 = p(6);
p_g2 = p(7);

p_c  = p(8);

% ------------------------------------------------------------
% 3. Inertia matrix M(q)
% ------------------------------------------------------------
c2 = cos(theta2);
s2 = sin(theta2);

M11 = p_a + 1 + 2*p_c*c2;
M12 = 1 + p_c*c2;
M21 = M12;
M22 = 1;

% ------------------------------------------------------------
% 4. Coriolis/centrifugal term C(q,qdot)*qdot
% ------------------------------------------------------------
Cqdot1 = -p_c*s2*(2*dtheta1*dtheta2 + dtheta2^2);
Cqdot2 =  p_c*s2*dtheta1^2;

% ------------------------------------------------------------
% 5. Viscous friction vector F(qdot)
% ------------------------------------------------------------
F1 = p_b1*dtheta1;
F2 = p_b2*dtheta2;

% ------------------------------------------------------------
% 6. Gravity vector G(q)
% ------------------------------------------------------------
G1 = -p_g1*sin(theta1) - p_g2*sin(theta1 + theta2);
G2 = -p_g2*sin(theta1 + theta2);

% ------------------------------------------------------------
% 7. Input vector I(u)
% ------------------------------------------------------------
I1 = p_u*u + p_0;
I2 = 0;

% ------------------------------------------------------------
% 8. Solve M(q)*qddot = I - Cqdot - F - G
% ------------------------------------------------------------
rhs1 = I1 - Cqdot1 - F1 - G1;
rhs2 = I2 - Cqdot2 - F2 - G2;

detM = M11*M22 - M12*M21;

ddtheta1 = ( M22*rhs1 - M12*rhs2) / detM;
ddtheta2 = (-M21*rhs1 + M11*rhs2) / detM;

% ------------------------------------------------------------
% 9. State derivative
% ------------------------------------------------------------
xdot = zeros(4,1);

xdot(1) = dtheta1;
xdot(2) = dtheta2;
xdot(3) = ddtheta1;
xdot(4) = ddtheta2;

% ------------------------------------------------------------
% 10. Outputs
% ------------------------------------------------------------
y = zeros(2,1);

y(1) = theta1;
y(2) = theta2;

end
