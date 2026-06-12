clear; clc; close all;

% ============================================================
% simulate_hinf.m
% ------------------------------------------------------------
% Closes the loop between the discrete H-inf controller (Kd)
% and the nonlinear pendulum plant.  No Simulink needed.
%
% Prerequisites:
%   hinf_controller.mat         (from hinf_design.m)
%   linearized_uncertain_plant.mat  (from linearize_uncertain_plant.m)
% ============================================================

%% 1. Load controller
load('hinf_controller.mat', 'Ak', 'Bk', 'Ck', 'Dk', 'h');
nK = size(Ak, 1);
fprintf('Controller: %d states,  Ts = %.3f s\n', nK, h);

%% 2. Load plant parameters
load('linearized_uncertain_plant.mat', 'x0', 'P');
nom     = cell2mat(P(:,2));                    % 14x1 nominal parameter vector
csmooth = [nom(12); nom(13); nom(11)];         % [epsv1; epsv2; vs2]  must stay numeric

%% 3. Simulation settings  <-- edit here
Tsim = 10;                     % total simulation time [s]

% Operating point and reference
x_op = x0;                     % [pi; pi; 0; 0]
r    = x_op(1:2);              % hold at operating point [pi; pi]

% Initial condition: perturbation away from the operating point.
% The controller should bring it back.  Increase to stress-test.
dx0    = [0.15; 0.10; 0; 0];   % [rad; rad; rad/s; rad/s]  ~8.6° / ~5.7°

% Optional theta1 reference step (tests tracking as well as balancing)
do_ref_step   = true;
t_step        = 4.0;           % time at which step is applied [s]
ref_step_size = 0.15;          % step size [rad]  (~8.6°)

%% 4. Pre-allocate storage
t      = (0:h:Tsim)';
N      = length(t);
x_hist = zeros(4, N);
u_hist = zeros(1, N);
r_hist = zeros(2, N);

%% 5. Simulate
x  = x_op + dx0;               % initial plant state
xK = zeros(nK, 1);             % controller state (starts at zero)

for k = 1:N
    % --- Reference (apply step if requested) ---
    r_now = r;
    if do_ref_step && t(k) >= t_step
        r_now(1) = r(1) + ref_step_size;
    end

    % --- Measurements: direct angle readings ---
    y = x(1:2);

    % --- Error: [r1 - theta1 ; r2 - theta2] ---
    e = r_now - y;

    % --- Controller output ---
    u = Ck * xK + Dk * e;

    % --- Clamp to hardware voltage limits [-1, +1] V ---
    u = max(-1, min(1, u));

    % --- Store ---
    x_hist(:, k) = x;
    u_hist(k)    = u;
    r_hist(:, k) = r_now;

    % --- Advance controller state (forward Euler on discrete model) ---
    xK = Ak * xK + Bk * e;

    % --- Integrate nonlinear plant over one sample period (RK4) ---
    x = rk4_step(x, u, nom, csmooth, h);
end

%% 6. Plot results
dev1    = x_hist(1,:) - x_op(1);  % theta1 deviation from OP  [rad]
dev2    = x_hist(2,:) - x_op(2);  % theta2 deviation from OP  [rad]
ref_dev = r_hist(1,:) - x_op(1);  % theta1 reference deviation [rad]
sat_frac = mean(abs(u_hist) > 0.99) * 100;

figure('Name','H-inf nonlinear simulation','Position',[80 80 950 750]);

subplot(3,1,1);
plot(t, rad2deg(dev1), 'b', 'LineWidth',1.2); hold on;
plot(t, rad2deg(ref_dev), 'r--', 'LineWidth',1.0);
yline(0, 'k:', 'LineWidth', 0.8);
ylabel('\theta_1 deviation (deg)');
legend('\theta_1 (actual)', 'reference', 'Location','northeast');
grid on;
title(sprintf(['H-\\infty nonlinear closed-loop  ' ...
    '(IC: \\theta_1{+}%.1f\\circ, \\theta_2{+}%.1f\\circ)'], ...
    rad2deg(dx0(1)), rad2deg(dx0(2))));

subplot(3,1,2);
plot(t, rad2deg(dev2), 'b', 'LineWidth',1.2); hold on;
yline(0, 'r--', 'LineWidth',1.0);
ylabel('\theta_2 deviation (deg)');
legend('\theta_2 (actual)', 'target 0', 'Location','northeast');
grid on;

subplot(3,1,3);
plot(t, u_hist, 'b', 'LineWidth',1.2); hold on;
yline( 1, 'r--', 'LineWidth',1.0);
yline(-1, 'r--', 'LineWidth',1.0);
ylabel('u (V)'); xlabel('Time (s)');
legend('u', '{\pm}1 V limits', 'Location','northeast');
grid on;
fprintf('Control saturation: %.1f%% of samples hit the ±1 V limit\n', sat_frac);

%% ============================================================
%%  What to look for in the plots
%  theta1, theta2:  Should return to 0 deviation after the transient.
%  u:               Should mostly stay within ±1 V; heavy saturation
%                   means the initial perturbation is too large for the
%                   linear controller to handle or the gains are too high.
%  Reference step:  With eps_S = 0.5 there will be visible SS error in
%                   theta1 (expected — reduce eps_S in hinf_design.m to fix).
%  Divergence:      If theta1 or theta2 grows without bound the linear
%                   controller is failing on the nonlinear plant; try
%                   reducing dx0 to confirm the controller works near the OP.
%% ============================================================

% ============================================================
% Local functions
% ============================================================

function x_new = rk4_step(x, u, p, csmooth, h)
% 4th-order Runge–Kutta integration of the nonlinear plant over one step.
    f1 = nl_plant(x,            u, p, csmooth);
    f2 = nl_plant(x + h/2*f1,  u, p, csmooth);
    f3 = nl_plant(x + h/2*f2,  u, p, csmooth);
    f4 = nl_plant(x + h*f3,    u, p, csmooth);
    x_new = x + (h/6) * (f1 + 2*f2 + 2*f3 + f4);
end

function xdot = nl_plant(x, u, p, csmooth)
% Nonlinear pendulum ODE.  p is the 14x1 numeric parameter vector.
% csmooth = [epsv1; epsv2; vs2] — must be numeric scalars.
    theta1  = x(1);  theta2  = x(2);
    dtheta1 = x(3);  dtheta2 = x(4);

    pa=p(1); pb1=p(2); pc1=p(3); pg1=p(4); pu=p(5); p0=p(6);
    pb2=p(7); pg2=p(8); pc2=p(9); psdelta2=p(10);
    pc=p(14);

    epsv1 = csmooth(1);
    epsv2 = csmooth(2);
    vs2   = csmooth(3);

    c2 = cos(theta2);  s2 = sin(theta2);

    M11 = pa + 1 + 2*pc*c2;
    M12 = 1  +   pc*c2;
    M22 = 1;

    Cq1 = -pc*s2*(2*dtheta1*dtheta2 + dtheta2^2);
    Cq2 =  pc*s2*dtheta1^2;

    F1 = pb1*dtheta1 + pc1*tanh(dtheta1/epsv1);
    F2 = pb2*dtheta2 + (pc2 + psdelta2*exp(-(dtheta2/vs2)^2))*tanh(dtheta2/epsv2);

    G1 = -pg1*sin(theta1) - pg2*sin(theta1+theta2);
    G2 = -pg2*sin(theta1+theta2);

    I1 = pu*u + p0;

    rhs1 = I1 - Cq1 - F1 - G1;
    rhs2 =    - Cq2 - F2 - G2;

    detM = M11*M22 - M12^2;
    xdot = [dtheta1;
            dtheta2;
            ( M22*rhs1 - M12*rhs2) / detM;
            (-M12*rhs1 + M11*rhs2) / detM];
end