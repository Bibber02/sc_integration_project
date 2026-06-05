% ------------------------------------------------------------------------
% Structured H-infinity (closed-loop shaping) design

%% 1. Load the linearized plant ------------------------------------------
load('linearized_plant.mat','sys_lin');     % -> sys_lin, A,B,C,D, etc.

G = sys_lin;
G.InputName  = 'u';
G.OutputName = {'y1','y2'};                  % y1 = theta1, y2 = theta2

% Quick sanity check on the open-loop plant
ol_poles = pole(G);
fprintf('Open-loop poles:\n'); disp(ol_poles);
if any(real(ol_poles) > 0)
    fprintf('--> Plant is UNSTABLE (expected for the upright equilibrium).\n');
    fprintf('    Fastest unstable pole: %.3f rad/s\n', ...
            max(real(ol_poles(real(ol_poles)>0))));
    fprintf('    Choose wS comfortably ABOVE this value.\n\n');
end

%% 2. Define the tunable controller structure ----------------------------
C1 = tunablePID('C1','pid');   % acts on e1 = r1 - theta1  (primary loop)
C2 = tunablePID('C2','pd');         % PD only (no redundant integrator)

C1.u = 'e1';  C1.y = 'uc1';
C2.u = 'e2';  C2.y = 'uc2';

%% 3. Build the closed-loop interconnection (S/KS standard form) ----------
% External inputs : references r1,r2 and output disturbances do1,do2
% External outputs: measured y1,y2 ; control u ; errors e1,e2
Sum_y1 = sumblk('y1m = y1 + do1');     % measured = plant out + disturbance
Sum_y2 = sumblk('y2m = y2 + do2');
Sum_e1 = sumblk('e1  = r1 - y1m');
Sum_e2 = sumblk('e2  = r2 - y2m');
Sum_u  = sumblk('u   = uc1 + uc2');

CL0 = connect(G, C1, C2, Sum_y1, Sum_y2, Sum_e1, Sum_e2, Sum_u, ...
              {'r1','r2','do1','do2'}, ...        % external inputs
              {'y1m','y2m','u','e1','e2'});       % external outputs
% In CL0:  do -> [y1m;y2m] is the output sensitivity  S_o (2x2)
%          do -> u          is the control sensitivity -K*S_o (1x2)

%% 4. Weighting filters (inverses of the desired closed-loop shapes) ------
% --- Sensitivity weight  W_S = 1/S_ref  (low gain at HF, high gain at LF) ---
As = db2mag(-40);    % allowed LF sensitivity  (-> small steady-state error)   [knob]
Ms = 1.5;            % allowed peak sensitivity (~ GM>=8.6dB, PM>=37deg)       [knob]
wS = 5;              % target disturbance-rejection bandwidth [rad/s]   [MAIN knob]
WS = makeweight(1/As, [wS 1], 1/Ms);   % crosses 0 dB at ~wS

% --- Control-sensitivity weight  W_KS = 1/KS_ref (high gain at HF) ---
MU  = db2mag(40);    % LF control authority allowed (light LF penalty = 1/MU)   [knob]
Au  = db2mag(-20);   % HF control attenuation requested (heavy HF penalty=1/Au) [knob]
wKS = 50;            % control roll-off frequency [rad/s] (~ actuator bandwidth)[knob]
WKS = makeweight(1/MU, [wKS 1], 1/Au);

%% 5. Tuning goals -------------------------------------------------------
% Soft goals (minimized): shape S and KS. systune handles them separately,
% so each returns its own performance level (~gamma) -> aim for ~1.
R_S  = TuningGoal.WeightedGain({'do1','do2'},{'y1m','y2m'}, WS , []);
R_S.Name  = 'Sensitivity (S)';

R_KS = TuningGoal.WeightedGain({'do1','do2'}, 'u',          WKS, []);
R_KS.Name = 'Control sensitivity (KS)';

% Hard goal: keep ALL closed-loop poles stable and in a sensible region.
% (This is the recommended way to stop systune pushing poles to crazy
%  frequencies; tune MaxFreq near your actuator bandwidth.)
minDecay   = 1e-3;          % > 0  -> guarantees stability with a margin
minDamping = 0.05;          % light minimum damping                       [knob]
maxFreq    = 10*wKS;        % cap fastest closed-loop pole                 [knob]
R_Poles = TuningGoal.Poles(minDecay, minDamping, maxFreq);
R_Poles.Name = 'CL pole region';

SoftGoals = [R_S, R_KS];
HardGoals =  R_Poles;

%% 6. Tune with systune --------------------------------------------------
opt = systuneOptions('RandomStart', 5, 'UseParallel', false); % more starts = more robust
rng(0);  % reproducibility
[CL, fSoft, gHard, info] = systune(CL0, SoftGoals, HardGoals, opt);

fprintf('\n==================== systune results ====================\n');
fprintf('Soft-goal performance levels (gamma, target ~1):\n');
fprintf('   S  goal : %.4f\n', fSoft(1));
fprintf('   KS goal : %.4f\n', fSoft(2));
fprintf('Hard-goal value (must be <= 1): %.4f\n', gHard);
fprintf('=========================================================\n\n');

% Tuned controller blocks
C1t = getBlockValue(CL,'C1');
C2t = getBlockValue(CL,'C2');
fprintf('Tuned C1 (theta1 loop):\n');  disp(tf(C1t));
fprintf('Tuned C2 (theta2 loop):\n');  disp(tf(C2t));
showTunable(CL);

%% 7. Reconstruct and validate the actual closed loop --------------------
% Full 1x2 controller: u = [C1 C2] * [e1; e2]
Kc = [C1t , C2t];                 % 1x2 system mapping [e1;e2] -> u
Kc.InputName  = {'e1','e2'};
Kc.OutputName = 'u';

Lo  = G*Kc;                       % output loop gain (2x2)
So  = feedback(eye(2), Lo);       % output sensitivity  = inv(I+Lo)
To  = eye(2) - So;                % output co-sensitivity
KSo = Kc*So;                      % control sensitivity (1x2)

% Closed-loop stability check
cl_poles = pole(So);
fprintf('Closed-loop poles (real parts should all be < 0):\n');
disp(cl_poles);
if all(real(cl_poles) < 0)
    fprintf('--> Closed loop is STABLE.\n\n');
else
    warning('Closed loop is NOT stable - revisit the weights/goals.');
end

% Achieved peaks
fprintf('Achieved ||S||_inf  = %.3f  (peak sensitivity / margins)\n', hinfnorm(So));
fprintf('Achieved ||T||_inf  = %.3f\n', hinfnorm(To));

%% 8. Frequency-domain plots: gang-of-four vs the weighting targets -------
w = logspace(-2, 3, 600);

figure('Name','Closed-loop shaping');
subplot(2,2,1)
sigma(So, 'b', inv(WS), 'r--', w); grid on
title('Sensitivity S_o vs 1/W_S (target)'); legend('S_o','1/W_S','Location','SouthEast');

subplot(2,2,2)
sigma(KSo, 'b', inv(WKS), 'r--', w); grid on
title('Control sensitivity KS_o vs 1/W_{KS}'); legend('KS_o','1/W_{KS}','Location','NorthEast');

subplot(2,2,3)
sigma(To, 'b', w); grid on
title('Co-sensitivity T_o');

subplot(2,2,4)
sigma(Lo, 'b', w); grid on
title('Output loop gain L_o = G K');

%% 9. Time-domain check (linear) -----------------------------------------
% Response to a unit step output disturbance on theta1 (regulation task).
T_do2y = So;                      % do -> y
figure('Name','Disturbance rejection (linear)');
step(T_do2y(:,1), 5); grid on
title('Linear response of [\theta_1;\theta_2] to a step disturbance on \theta_1');

%% 10. Save the design ---------------------------------------------------
save('hinf_pid_controller.mat','Kc','C1t','C2t','So','To','KSo','Lo', ...
     'WS','WKS','fSoft','gHard','G');
fprintf('Saved controller and analysis to hinf_pid_controller.mat\n');
