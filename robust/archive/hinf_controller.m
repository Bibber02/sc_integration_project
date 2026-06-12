load('linearized_uncertain_plant.mat','sys_lin');

G = ss(sys_lin.NominalValue);
G.InputName  = 'u';
G.OutputName = {'theta1','theta2'};

ny = 2;   % outputs
nu = 1;   % inputs

%% Weighting filters
% W_S
% theta1 channel
M1   = 2.0;     
wb1 = 2.0;    
eps1 = 1e-3;
% theta2 channel
M2   = 2.0;     
wb2 = 4.0;    
eps2 = 1e-4;

Ws1 = tf([1/M1  wb1],[1  wb1*eps1]);
Ws2 = tf([1/M2  wb2],[1  wb2*eps2]);
WS  = blkdiag(Ws1, Ws2);

% W_KS
Mu   = 30;      
wbc = 50;     
epsu = 1e-3;

Wu   = tf([1  wbc/Mu],[epsu  wbc]);
WKS  = Wu;                                   

% W_T
Mt   = 2.0;     
wbt = 30;     
epst = 1e-3;

Wt   = tf([1  wbt/Mt],[epst*1  wbt]); 
WT   = blkdiag(Wt, Wt);

%% Generalized plant
P = augw(G, WS, WKS);          
%P = augw(G, WS, WKS, WT);

%% Hinf Synthesis
[K, CL, gamma] = hinfsyn(P, ny, nu, 'Method','ric','Display','on');
fprintf('Achieved gamma = %.4f\n', gamma);
if gamma > 1.2
    warning(['gamma = %.2f is high'], gamma);
end

%% Closed-loop and checks
K.InputName  = {'theta1','theta2'};
K.OutputName = 'u';

L  = G*K;                       
So = feedback(eye(ny), G*K); 
To = eye(ny) - So;

figure; sigma(So, 'b', To, 'r', {1e-1,1e3}); grid on;
legend('S','T'); title('Output sensitivity and complementary sensitivity');

% Damping check
disp('Closed-loop poles:'); damp(feedback(G,K));

% plotz
figure; 
impulse(feedback(K, G), 5);
title('Controller output to disturbance (watch |u| <= 1)');

figure; 
step(feedback(K, G), 5);
title('Controller output to step');