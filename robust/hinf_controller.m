clear; clc; close all;

% ============================================================
% hinf_design.m
% ------------------------------------------------------------
% Designs a full-order H-inf output-feedback controller for
% the rotational pendulum, then discretizes it for use as an
% LTI System block in Simulink.
%
% Prerequisites:
%   - Run linearize_uncertain_plant.m first to generate
%     linearized_uncertain_plant.mat
% ============================================================

%% 1. Load nominal plant
load('linearized_uncertain_plant.mat', 'sys_lin');
Gnom = sys_lin.NominalValue;        % ss, 2 outputs (theta1,theta2), 1 input (u)

fprintf('Nominal open-loop eigenvalues:\n');
ev = eig(Gnom.A);
disp(ev);
fprintf('Unstable poles: %d  (your bandwidth must exceed their magnitude)\n\n', ...
        sum(real(ev) > 0));

% ============================================================
%% 2. TUNING PARAMETERS  <-- edit here
% ============================================================
h        = 0.01;    % hardware sample period [s] — match your Simulink model
omega_B  = 8;      % target closed-loop bandwidth [rad/s]
                    % rule of thumb: at least 3-5x the unstable pole magnitude
M_S      = 2.0;     % max sensitivity peak |S|_inf (1.2 to 2 is typical)
eps_S    = 0.1;    % near-integrator factor; |W_S(0)| ~ omega_B/(eps_S*M_S)
                    % smaller = tighter DC regulation, but harder to achieve
omega_KS = 200;     % control effort rolloff [rad/s], keep below h^-1/2
M_KS     = 5.0;     % high-frequency bound on |KS|; increase if synthesis fails

% Independent tuning per output channel (set equal to startxsxxsssssxxs, then split)
omega_B2 = omega_B; % bandwidth weight for theta2 channel
M_S2     = M_S;
% ============================================================

%% 3. Weighting filters
s = tf('s');

% W_S(s) = (s/M_S + omega_B) / (s + omega_B*eps_S)
%   |W_S(jw)| ~ 1/M_S  at HF  ->  |S| can peak up to M_S
%   |W_S(0)|  = omega_B/(eps_S*M_S) >> 1  ->  tight DC regulation
%   Crossover near omega_B  ->  sets the bandwidth target
W_s1 = (s/M_S  + omega_B ) / (s + omega_B *eps_S);   % theta1 channel
W_s2 = (s/M_S2 + omega_B2) / (s + omega_B2*eps_S);   % theta2 channel
W_S  = blkdiag(W_s1, W_s2);                           % 2x2 diagonal

% W_KS(s): roll off controller gain beyond actuator bandwidth
%   Low gain at LF -> allow control effort; high gain at HF -> limit it
%   u is bounded to [-1,+1] V on the hardware
W_KS = (s + omega_KS/M_KS) / (M_KS * (s + omega_KS));

%% 4. Augmented plant and H-inf synthesis
% augw(G, W1, W2, []) sets up the standard S/KS mixed-sensitivity problem:
%   minimise  ||[W_S * S_o; W_KS * K*S_o]||_inf  over stabilising K
P = augw(Gnom, W_S, W_KS, []);     % nmeas = 2 (outputs), ncont = 1 (input)

%opts  = hinfsynOptions('Method', 'LMI');    % or 'RIC' for Riccati-based
[K, ~, gamma] = hinfsyn(P, 2, 1);

fprintf('H-inf synthesis complete.\n');
fprintf('  gamma achieved = %.4f\n', gamma);
if gamma < 1
    fprintf('  gamma < 1: all weighted targets met.\n\n');
else
    warning(['gamma > 1: weighted targets not fully met.\n' ...
             'Try: reducing omega_B, increasing M_S/M_KS, or ' ...
             'increasing eps_S.']);
end

%% 5. Nominal closed-loop analysis
F = loopsens(Gnom, K);
fprintf('Nominal closed-loop stable: %d  (1 = yes)\n\n', F.Stable);
if ~F.Stable
    error('Nominal closed loop is UNSTABLE. Revisit the weights.');
end

% Singular value plot of nominal S and T
figure('Name','Nominal sensitivity');
sigma(F.So, 'b', F.To, 'r', {1e-1, 1e3});
legend('S_o','T_o'); grid on;
title(sprintf('Nominal closed-loop singular values  (\\gamma = %.3f)', gamma));
xlabel('Frequency (rad/s)'); ylabel('Singular values (dB)');

% Closed-loop step response (reference tracking)
figure('Name','Nominal step response');
step(F.To, 5); grid on;
title('Nominal closed-loop step response (T_o)');

%% 6. Robustness check on the UNCERTAIN plant
% K was designed on the nominal; verify it survives the full uncertainty set.
% This is the payoff for building sys_lin as a uss earlier.
fprintf('--- Robust stability analysis ---\n');
Func = loopsens(sys_lin, K);        % sys_lin is the uss
opt  = robOptions('Display', 'on', 'Sensitivity', 'on');
[SM, wcu] = robstab(Func.So, opt);
fprintf('Robust stability margin: LB = %.3f,  UB = %.3f\n', ...
        SM.LowerBound, SM.UpperBound);
fprintf('(margin > 1 -> stable for all modelled uncertainty)\n\n');

% Worst-case gain of output sensitivity (indicator of worst-case overshoot)
[wcg, ~] = wcgain(Func.So);
fprintf('Worst-case sensitivity peak: LB = %.3f,  UB = %.3f\n\n', ...
        wcg.LowerBound, wcg.UpperBound);

%% 7. Discretise for Simulink
% Tustin (bilinear) preserves phase behaviour near crossover and is
% standard for output-feedback controllers.
Kd = c2d(K, h, 'tustin');
Kd.InputName  = {'e_theta1', 'e_theta2'};
Kd.OutputName = {'u'};

fprintf('Discrete controller: %d states,  Ts = %.4f s\n\n', order(Kd), h);

% Individual matrices — use these if you prefer a Discrete State-Space block
Ak = Kd.A;  Bk = Kd.B;  Ck = Kd.C;  Dk = Kd.D;

%% 8. Save
save('hinf_controller.mat', ...
     'Kd', 'Ak', 'Bk', 'Ck', 'Dk', ...
     'K', 'gamma', 'h', 'omega_B', 'omega_B2', 'W_S', 'W_KS');

disp('=== DONE ===');
disp('Controller saved to hinf_controller.mat.');
disp(' ');
disp('Simulink (LTI System block):');
disp('  1. Library Browser -> Control System Toolbox -> LTI System');
disp('  2. Block parameter "LTI system variable" = Kd');
disp('  3. Input: 2x1 Mux  ->  [r1 - theta1 ;  0 - theta2]');
disp('  4. Output: u -> motor voltage');
disp(' ');
disp('Simulink (Discrete State-Space block, no toolbox block needed):');
disp('  A=Ak  B=Bk  C=Ck  D=Dk  Ts=h  (all in workspace after this script)');