clear;
clc;
close all;

% ------------------------------------------------------------
% 1. Load parameter vector p
% ------------------------------------------------------------
% This should create p in the workspace.
% p = [pa; pb1; pc1; pg1; pu; p0; pb2; pg2; pc2; psdelta2; vs2; epsv1; epsv2; pc]
load_parameters;

% ------------------------------------------------------------
% 2. Define operating point
% ------------------------------------------------------------
% State order:
% x = [theta1; theta2; theta1_dot; theta2_dot]

x0 = [pi; pi; 0; 0];

% Input equation:
% I1 = pu*u + p0
%
% Therefore, to cancel the torque offset at equilibrium:
% pu*u0 + p0 = 0
%
% p(5) = pu
% p(6) = p0

u0 = -p(6)/p(5);

% ------------------------------------------------------------
% 3. Check nonlinear model at operating point
% ------------------------------------------------------------

[f0, y0] = nonlinearPlant(x0, u0, p);

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

% Finite-difference step sizes
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

    [f_plus,  y_plus]  = nonlinearPlant(x0 + dx, u0, p);
    [f_minus, y_minus] = nonlinearPlant(x0 - dx, u0, p);

    A(:,i) = (f_plus - f_minus)/(2*hx);
    C(:,i) = (y_plus - y_minus)/(2*hx);
end

% ------------------------------------------------------------
% 6. Compute B and D matrices
% ------------------------------------------------------------
% B = df/du
% D = dy/du

[f_plus,  y_plus]  = nonlinearPlant(x0, u0 + hu, p);
[f_minus, y_minus] = nonlinearPlant(x0, u0 - hu, p);

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

Ts = 0.001;   % Sampling period [s]

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

disp('Continuous-time model:');
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
disp('Discrete-time model:');
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

save('linearized_plant.mat', ...
     'A','B','C','D', ...
     'Ad','Bd','Cd','Dd', ...
     'sys_lin','sys_disc', ...
     'Ts','x0','u0','p');

disp('Saved continuous and discrete linearized models to linearized_plant.mat');


function [xdot, y] = nonlinearPlant(x, u, p)
%#codegen
% nonlinearPlant

% Reduced nonlinear plant model for the rotational double pendulum.
%
% Inputs:
%   x = [theta1; theta2; theta1_dot; theta2_dot]
%   u = motor input command / voltage-like input
%   p = parameter vector loaded from .mat file
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
pa       = p(1);

pb1      = p(2);
pc1      = p(3);
pg1      = p(4);
pu       = p(5);
p0       = p(6);

pb2      = p(7);
pg2      = p(8);
pc2      = p(9);
psdelta2 = p(10);
vs2      = p(11);

epsv1    = p(12);
epsv2    = p(13);

pc       = p(14);

% ------------------------------------------------------------
% 3. Inertia matrix M(q)
% ------------------------------------------------------------
c2 = cos(theta2);
s2 = sin(theta2);

M11 = pa + 1 + 2*pc*c2;
M12 = 1 + pc*c2;
M21 = M12;
M22 = 1;

% ------------------------------------------------------------
% 4. Coriolis/centrifugal term C(q,qdot)*qdot
% ------------------------------------------------------------
Cqdot1 = -pc*s2*(2*dtheta1*dtheta2 + dtheta2^2);
Cqdot2 =  pc*s2*dtheta1^2;

% ------------------------------------------------------------
% 5. Friction vector F(qdot)
% ------------------------------------------------------------
F1 = pb1*dtheta1 + pc1*tanh(dtheta1/epsv1);

stribeck2 = pc2 + psdelta2*exp(-(dtheta2/vs2)^2);
F2 = pb2*dtheta2 + stribeck2*tanh(dtheta2/epsv2);

% ------------------------------------------------------------
% 6. Gravity vector G(q)
% ------------------------------------------------------------
G1 = -pg1*sin(theta1) - pg2*sin(theta1 + theta2);
G2 = -pg2*sin(theta1 + theta2);

% ------------------------------------------------------------
% 7. Input vector I(u)
% ------------------------------------------------------------
I1 = pu * u + p0;
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