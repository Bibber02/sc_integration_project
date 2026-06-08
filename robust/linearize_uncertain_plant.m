clear;
clc;
close all;

% ============================================================
% linearize_uncertain_plant.m
% ============================================================

% ------------------------------------------------------------
% Parameter table  --  EDIT THESE
% ------------------------------------------------------------
% Order MUST match the unpacking inside nonlinearPlant:
%   [pa; pb1; pc1; pg1; pu; p0; pb2; pg2; pc2; psdelta2; vs2; epsv1; epsv2; pc]
% Columns: { name , nominal value , standard deviation sigma }
%   sigma = 0  ->  treated as CERTAIN (fixed, no ureal created)

P = {
%   name        nominal       sigma
    'pa',        64.7915,      0.5849
    'pb1',       3.2117e+03,   3.0552
    'pc1',       635.3865,     2.0833
    'pg1',       430.5226,     0.5236
    'pu',        2.2424e+04,   14.3848
    'p0',        0,     0    % NOTE: drops out of A,B,C,D (constant offset)
    'pb2',       0.0670,       0          % passive link: locked
    'pg2',       112.0268,     0         % locked (add sigma if you have one)
    'pc2',       0.2491,       0          % locked
    'psdelta2',  0.1836,       0         % locked
    'vs2',       3.5000,       0          % MUST stay fixed (enters exp); also inert at this OP
    'epsv1',     0.0500,       0          % MUST stay fixed (enters tanh)
    'epsv2',     0.0300,       0          % MUST stay fixed (enters tanh)
    'pc',        0.5300,       0.05          % <-- PLACEHOLDER. pc = (l1/g)*pg2.
};

% ------------------------------------------------------------
% Uncertainty settings
% ------------------------------------------------------------
kSigma  = 2;                               % uncertainty range = nominal +/- kSigma*sigma
mustFix = {'vs2','epsv1','epsv2'};         % never make these uncertain (inside tanh/exp)

% ------------------------------------------------------------
% Build the parameter vector p
% ------------------------------------------------------------
names = P(:,1);
nom   = cell2mat(P(:,2));
sig   = cell2mat(P(:,3));
np    = numel(names);

plist    = cell(np,1);
uncNames = {};
for i = 1:np
    if sig(i) > 0 && ~ismember(names{i}, mustFix)
        plist{i} = ureal(names{i}, nom(i), 'Plusminus', kSigma*sig(i));
        uncNames{end+1} = names{i}; %#ok<SAGROW>
    else
        plist{i} = nom(i);
    end
end
p = vertcat(plist{:});      % becomes a umat column if any ureal is present

fprintf('Uncertain parameters (%d): %s\n', numel(uncNames), strjoin(uncNames, ', '));

% ------------------------------------------------------------
% Fixed smoothing / Stribeck constants -> kept strictly NUMERIC
% -----------------------f-------------------------------------
% These are the constants that appear inside tanh()/exp(). Indexing a
% umat returns a umat even for constant entries, and tanh/exp are not
% defined for umat. So we extract them here as plain doubles from the
% NOMINAL table and pass them to nonlinearPlant separately.
%   csmooth = [epsv1; epsv2; vs2]
csmooth = [nom(12); nom(13); nom(11)];

% ------------------------------------------------------------
% Operating point (computed from NOMINAL values -> numeric)
% ------------------------------------------------------------
% State order: x = [theta1; theta2; theta1_dot; theta2_dot]
x0 = [pi; pi; 0; 0];

% Equilibrium input cancels the torque offset: pu*u0 + p0 = 0
u0 = -nom(6)/nom(5);        % nominal -> numeric, so the OP is fixed

% Equilibrium check on the NOMINAL model
[f0nom, y0] = nonlinearPlant(x0, u0, nom, csmooth);
fprintf('||f(x0,u0)||_nominal = %.3e\n', norm(f0nom));
if norm(f0nom) > 1e-6
    warning(['The operating point is not an exact equilibrium for the ', ...
             'nominal parameters; the linear model has a small affine offset.']);
end

% ------------------------------------------------------------
% Numerical linearization settings
% ------------------------------------------------------------
nx = numel(x0);
nu = 1;
ny = numel(y0);

hx = 1e-6;
hu = 1e-6;

% ------------------------------------------------------------
% A = df/dx and C = dy/dx
% ------------------------------------------------------------
Acols = cell(1, nx);
Ccols = cell(1, nx);
for i = 1:nx
    dx = zeros(nx,1);
    dx(i) = hx;

    [f_plus,  y_plus ] = nonlinearPlant(x0 + dx, u0, p, csmooth);
    [f_minus, y_minus] = nonlinearPlant(x0 - dx, u0, p, csmooth);

    Acols{i} = (f_plus  - f_minus )/(2*hx);
    Ccols{i} = (y_plus  - y_minus )/(2*hx);
end
A = [Acols{:}];
C = [Ccols{:}];

% ------------------------------------------------------------
% B = df/du and D = dy/du
% ------------------------------------------------------------
[f_plus,  y_plus ] = nonlinearPlant(x0, u0 + hu, p, csmooth);
[f_minus, y_minus] = nonlinearPlant(x0, u0 - hu, p, csmooth);

B = (f_plus - f_minus)/(2*hu);
D = (y_plus - y_minus)/(2*hu);

% A = simplify(A,'full'); B = simplify(B,'full');
% C = simplify(C,'full'); D = simplify(D,'full');

% ------------------------------------------------------------
% Assemble the (uncertain) continuous-time state-space model
% ------------------------------------------------------------
sys_lin = ss(A, B, C, D);
sys_lin.StateName  = {'theta1','theta2','theta1_dot','theta2_dot'};
sys_lin.InputName  = {'u'};
sys_lin.OutputName = {'theta1','theta2'};

isUncertain = isa(sys_lin,'uss');

% ------------------------------------------------------------
% Discretize (ZOH)
% ------------------------------------------------------------
Ts = 0.001;
try
    sys_disc = c2d(sys_lin, Ts, 'zoh');
catch ME
    warning('c2d on the uss failed (%s). Discretizing the NOMINAL model instead.', ME.message);
    sys_disc = c2d(sys_lin.NominalValue, Ts, 'zoh');
end

% ------------------------------------------------------------
% Report
% ------------------------------------------------------------
disp(' ');
disp('Continuous-time uncertain model:');
sys_lin %#ok<NOPTS>

if isUncertain
    disp('Nominal continuous-time eigenvalues:');
    disp(eig(sys_lin.NominalValue.A));
else
    disp('Model is NOT uncertain (all sigma = 0). Eigenvalues:');
    disp(eig(A));
end

% ------------------------------------------------------------
% Visualization
% ------------------------------------------------------------
doPlots = true;
if doPlots && isUncertain
    N = 40;                               % number of random samples
    figure;
    subplot(1,1,1);
    bodemag(usample(sys_lin, N), 'b', sys_lin.NominalValue, 'r', {1e-1, 1e3});
    title('Bode magnitude: random samples (blue) vs nominal (red)');
    legend('samples','nominal');
    figure;
    ylim([-0.4 0.8]);
    legend('samples','nominal');
    subplot(1,1,1);
    step(usample(sys_lin, N), 'b', sys_lin.NominalValue, 'r', 10);
    title('Step response: random samples (blue) vs nominal (red)');
end

% ------------------------------------------------------------
% Save
% ------------------------------------------------------------
save('linearized_uncertain_plant.mat', ...
     'sys_lin','sys_disc','Ts','x0','u0','P','kSigma','uncNames');
disp('Saved uncertain linear model to linearized_uncertain_plant.mat');

% ============================================================
% Local function: nonlinear plant.
% Signature changed: the fixed smoothing/Stribeck constants are now
% passed as a separate NUMERIC vector csmooth = [epsv1; epsv2; vs2],
% so tanh()/exp() always receive double arguments even when p is a umat.
% ============================================================
function [xdot, y] = nonlinearPlant(x, u, p, csmooth)
%#codegen
theta1  = x(1);
theta2  = x(2);
dtheta1 = x(3);
dtheta2 = x(4);

% Uncertain / identified parameters (may be umat entries)
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
pc       = p(14);

% Fixed smoothing / Stribeck constants -> strictly numeric
epsv1    = csmooth(1);
epsv2    = csmooth(2);
vs2      = csmooth(3);

% Inertia matrix M(q)
c2 = cos(theta2);
s2 = sin(theta2);

M11 = pa + 1 + 2*pc*c2;
M12 = 1 + pc*c2;
M21 = M12;
M22 = 1;

% Coriolis/centrifugal C(q,qdot)*qdot
Cqdot1 = -pc*s2*(2*dtheta1*dtheta2 + dtheta2^2);
Cqdot2 =  pc*s2*dtheta1^2;

% Friction F(qdot)   (tanh/exp arguments are numeric -> OK with umat params)
F1 = pb1*dtheta1 + pc1*tanh(dtheta1/epsv1);

stribeck2 = pc2 + psdelta2*exp(-(dtheta2/vs2)^2);
F2 = pb2*dtheta2 + stribeck2*tanh(dtheta2/epsv2);

% Gravity G(q)
G1 = -pg1*sin(theta1) - pg2*sin(theta1 + theta2);
G2 = -pg2*sin(theta1 + theta2);

% Input I(u)
I1 = pu * u + p0;
I2 = 0;

% Solve M(q)*qddot = I - Cqdot - F - G
rhs1 = I1 - Cqdot1 - F1 - G1;
rhs2 = I2 - Cqdot2 - F2 - G2;

detM = M11*M22 - M12*M21;

ddtheta1 = ( M22*rhs1 - M12*rhs2) / detM;
ddtheta2 = (-M21*rhs1 + M11*rhs2) / detM;

xdot    = [dtheta1; dtheta2; ddtheta1; ddtheta2];
y       = [theta1; theta2];
end