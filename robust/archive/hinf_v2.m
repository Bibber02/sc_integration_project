%% Hinf REGULATION controller - rotational double pendulum
%  Underactuated: 1 input (motor on link 1), 2 measured outputs (theta1, theta2).
%  Setpoint: both links hanging down  [theta1 theta2 theta1dot theta2dot] = [pi 0 0 0].
%  Goal: REGULATE the equilibrium (reject input disturbance / initial-condition
%        energy), NOT track a reference. r = 0 throughout.
%
%  Output feedback: controller sees only the two angles (theta1, theta2),
%  same measurements as the LQG controller.
%
%  Key difference vs a fully-actuated 2x2 mixed-sensitivity setup:
%  we do NOT demand independent sensitivity shaping on both outputs (impossible
%  with one actuator). Instead we minimise ONE weighted regulated-error vector
%  against an input disturbance, and trade off theta1 vs theta2 error explicitly
%  via rho1/rho2 (the analogue of the LQG Q diagonal).

clear; clc;

%% ----------------------------------------------------------------------
%  Plant (linearised about the both-down equilibrium)
%  NOTE: load the plant linearised about [pi 0 0 0]. If your .mat was
%  linearised about a different operating point, re-linearise there first.
% ----------------------------------------------------------------------
load('linearized_uncertain_plant.mat','sys_lin');
G = ss(sys_lin.NominalValue);
G.InputName  = 'u';
G.OutputName = {'theta1','theta2'};

ny = 2;   % measured outputs fed to controller (theta1, theta2)
nu = 1;   % control inputs (motor)

% Sanity checks on the plant
fprintf('Plant: %d states, %d inputs, %d outputs\n', order(G), size(G,2), size(G,1));
ev = eig(G.A);
fprintf('Open-loop poles:\n'); disp(ev);
if any(real(ev) > 1e-6)
    warning('Plant has RHP poles - expected for "up" equilibria, NOT for both-down.');
end

%% ----------------------------------------------------------------------
%  Weighting filters  (course convention: makeweight(dcgain,[freq,mag],hfgain))
% ----------------------------------------------------------------------

% --- Performance weight We on the regulated error e = -y (since r=0) -----
%  Low-pass-inverse shape: high gain at LF (force error -> 0 in steady state),
%  rolls off at HF. 1/A_e is the demanded LF disturbance attenuation.
%  Relaxed vs the original (40 dB LF, not 80 dB) so a single actuator can meet it.
A_e   = 1e-2;     % LF gain of We ~ 1/A_e = 40 dB demanded attenuation
wb_e  = 3.0;      % desired closed-loop error bandwidth [rad/s]  <-- main knob
M_e   = 2.0;      % allowed sensitivity peak (Ms ~ 2 -> ~ GM 6 dB, PM ~ 29 deg)
We_scalar = makeweight(1/A_e, [wb_e, 1], 1/M_e);   % 1/We is the target Sref

% Per-angle trade-off (you CANNOT regulate both independently; pick the balance).
% Larger rho_i => tighter regulation demanded on that angle.
rho1 = 1.0;       % weight on theta1 error
rho2 = 1.0;       % weight on theta2 error  <-- raise if link 2 is the priority
We = blkdiag(rho1*We_scalar, rho2*We_scalar);
We.u = {'e1','e2'};  We.y = {'ze1','ze2'};

% --- Control-effort weight Wu on u  (proper high-pass, fixes the old Wu) --
%  High-pass: cheap control at LF, penalise control action beyond the
%  actuator bandwidth so |u| stays sane and HF actuator excitation is limited.
A_u   = 1e-2;     % LF gain (cheap control at low freq)
wb_u  = 50;       % actuator bandwidth / control roll-up freq [rad/s]
M_u   = 30;       % HF gain ceiling on the control penalty
Wu = makeweight(A_u, [wb_u, 1], M_u);
Wu.u = 'u';  Wu.y = 'zu';

%% ----------------------------------------------------------------------
%  Build the generalised plant P by interconnection (regulation, r = 0)
%
%   exogenous input :  di   (input disturbance, added to u at the plant input)
%   control input   :  u    (from controller)
%   regulated out   :  z  = [ ze1 ; ze2 ; zu ]   (We*e  and  Wu*u)
%   measured out    :  v  = [ e1 ; e2 ] = -[theta1 ; theta2]   (since r = 0)
%
%  Plant is driven by (u + di). Controller sees the negated angles (error to 0).
% ----------------------------------------------------------------------
sumIn = sumblk('up = u + di');          % up = plant input = control + disturbance
G.u = 'up';                              % plant input renamed

% measured error e = r - y = -y  (r = 0)
e1blk = sumblk('e1 = -theta1');
e2blk = sumblk('e2 = -theta2');
P = connect(G, sumIn, e1blk, e2blk, We, Wu, ...
            {'di','u'}, {'ze1','ze2','zu','e1','e2'});

nz = 3;   % regulated outputs (ze1, ze2, zu)
nv = 2;   % measurements to controller (e1, e2)

%% ----------------------------------------------------------------------
%  Hinf synthesis
% ----------------------------------------------------------------------
opts = hinfsynOptions('Method','lmi','Display','on');   % LMI is robust for hand-built P
[K, CL, gamma] = hinfsyn(P, nv, nu, opts);

fprintf('\nAchieved gamma = %.4f\n', gamma);
if gamma > 1.3
    warning('gamma = %.2f is high. Relax We (raise A_e / lower wb_e) or rebalance rho1/rho2.', gamma);
elseif gamma < 0.7
    fprintf('gamma < 0.7: spec is loose, you can tighten We (lower A_e / raise wb_e).\n');
else
    fprintf('gamma is in a healthy range (~1). Good starting design.\n');
end

%% ----------------------------------------------------------------------
%  Close the loop manually for analysis (note the sign: K acts on e = -y)
% ----------------------------------------------------------------------
K.u = {'e1','e2'};   K.y = 'u';

% Reconstruct the true feedback loop: u = K*(-y)
Kfb = K * [-1 0; 0 -1];          % maps y -> u directly (folds in the r=0 error sign)
Kfb = ss(Kfb);
Kfb.u = {'theta1','theta2'};  Kfb.y = 'u';

Graw = ss(sys_lin.NominalValue);
Graw.u = 'u';  Graw.y = {'theta1','theta2'};

L  = Graw*Kfb;                    % output loop gain Lo
So = feedback(eye(ny), L);        % output sensitivity
To = eye(ny) - So;                % complementary sensitivity
KS = feedback(Kfb, Graw);         % control sensitivity K*So

%% ----------------------------------------------------------------------
%  Checks / plots
% ----------------------------------------------------------------------
figure; sigma(So,'b', To,'r', {1e-2,1e3}); grid on;
legend('S_o','T_o'); title('Output sensitivity and complementary sensitivity');

figure; sigma(KS,'m', {1e-2,1e3}); grid on;
title('Control sensitivity K S_o  (watch HF roll-off)');

fprintf('\nClosed-loop poles (should all be in LHP):\n');
damp(feedback(Graw, Kfb));

% Disturbance-rejection time response: input-disturbance impulse -> angles
Tdy = feedback(Graw, Kfb);                     % di acts like an input on Graw path
figure; impulse(So, 5); grid on;
title('Output sensitivity impulse response (disturbance -> angles)');

% Control effort to a unit input disturbance (keep |u| within actuator limit)
figure; impulse(KS, 5); grid on;
title('Control effort to input disturbance  (watch |u|)');