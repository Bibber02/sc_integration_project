clear;
clear functions;
clc;
close all;

% --- 0. Start Console Logging ---
diary('estimation_console_output.txt');
fprintf('\n======================================================\n');
fprintf('Starting Full System Identification run...\n');
fprintf('======================================================\n\n');

% --- 1. Load PRBS Measurement Data ---
data_train = load("identification_data\prbs_train_121s_amp_015.mat");
data_val = load("identification_data\prbs_validation_121s_amp_015.mat");

Ts = 0.01; % Sampling period

t_train = data_train.theta_1.Time(:);
theta_1_train = data_train.theta_1.Data(:);
theta_2_train = data_train.theta_2.Data(:);
u_train = data_train.u_ts.Data(1:12001);

t_val = data_val.theta_1.Time(:);
theta_1_val = data_val.theta_1.Data(:);
theta_2_val = data_val.theta_2.Data(:);
u_val = data_val.u_ts.Data(1:12001);

% Make sure all signals have equal length
N_train = min([length(t_train), length(theta_1_train), length(theta_2_train), length(u_train)]);
N_val   = min([length(t_val),   length(theta_1_val),   length(theta_2_val),   length(u_val)]);

t_train       = t_train(1:N_train);
theta_1_train = theta_1_train(1:N_train);
theta_2_train = theta_2_train(1:N_train);
u_train       = -u_train(1:N_train);

t_val       = t_val(1:N_val);
theta_1_val = theta_1_val(1:N_val);
theta_2_val = theta_2_val(1:N_val);
u_val       = -u_val(1:N_val);

% --- 2. Implement Identified Passive Linkage Parameters ---
p_b2_val      = 0.067028846;
p_g2_val      = 112.02675;
p_c2_val      = 0.24908901;
p_sdelta2_val = 0.18360238;
v_s2_val      = 3.5;
eps_v2_val    = 0.03;

% --- 3. Format iddata & Truncate Validation Data ---
y_train = [theta_1_train, theta_2_train];
y_val = [theta_1_val, theta_2_val];

iddata_train = iddata(y_train, u_train, Ts);
iddata_val = iddata(y_val, u_val, Ts);

iddata_train.InputName  = {'u'}; iddata_train.OutputName = {'theta_1', 'theta_2'};
iddata_train.InputUnit  = {'V'}; iddata_train.OutputUnit = {'rad', 'rad'};
iddata_train.TimeUnit   = 's';

iddata_val.InputName  = {'u'};   iddata_val.OutputName = {'theta_1', 'theta_2'};
iddata_val.InputUnit  = {'V'};   iddata_val.OutputUnit = {'rad', 'rad'};
iddata_val.TimeUnit   = 's';

% Slice the validation data from t = 60s (index 6000) to the end
start_idx = 6000;
iddata_val = iddata_val(start_idx:end);

initial_states_train = {theta_1_train(1); theta_2_train(1); 0; 0};
initial_states_val   = {theta_1_val(start_idx); theta_2_val(start_idx); 0; 0};

% --- 4. Setup Grey-Box Model ---
model_file = 'greybox_id_full_model';
order = [2 1 4];     % 2 outputs, 1 input, 4 states
Ts_model = 0;        % continuous-time model

% 13 Parameters total (6 Active, 7 Passive/Smoothing)
parameters = {
    14;            % 1: p_a
    670;           % 2: p_b1
    10;            % 3: p_c1 (Active Coulomb)
    97;            % 4: p_g1
    4000;          % 5: p_u
    0;             % 6: p_0 (Torque Offset)
    p_b2_val;      % 7: p_b2 (Locked)
    p_g2_val;      % 8: p_g2 (Locked)
    p_c2_val;      % 9: p_c2 (Locked)
    p_sdelta2_val; % 10: p_sdelta2 (Locked)
    v_s2_val;      % 11: v_s2 (Locked)
    0.05;          % 12: eps_v1 (Active smoothing, Locked)
    eps_v2_val     % 13: eps_v2 (Locked)
};

model_train_initial = idnlgrey(model_file, order, parameters, initial_states_train, Ts_model);

parameter_names = {'p_a', 'p_b1', 'p_c1', 'p_g1', 'p_u', 'p_0', 'p_b2', 'p_g2', 'p_c2', 'p_sdelta2', 'v_s2', 'eps_v1', 'eps_v2'};

% [UPDATE]: p_c1 maximum bound increased to 250
minimum_values = [5, 100,  0,   5,  500, -200, 0, 0, 0, 0, 0, 0.001, 0.001];
maximum_values = [80, 1500, 250, 200, 12000, 200, Inf, Inf, Inf, Inf, Inf, 1, 1];

% Lock indices 7 through 13
fixed_flags = [false, false, false, false, false, false, true, true, true, true, true, true, true];

for k = 1:length(parameter_names)
    model_train_initial.Parameters(k).Name = parameter_names{k};
    model_train_initial.Parameters(k).Minimum = minimum_values(k);
    model_train_initial.Parameters(k).Maximum = maximum_values(k);
    model_train_initial.Parameters(k).Fixed = fixed_flags(k);
end

% Configure input/output/state names
model_train_initial.InputName = {'u'};  model_train_initial.InputUnit = {'V'};
model_train_initial.OutputName = {'theta_1', 'theta_2'}; model_train_initial.OutputUnit = {'rad', 'rad'};
model_train_initial.TimeUnit = 's';

model_train_initial.InitialStates(1).Name = 'theta_1'; model_train_initial.InitialStates(1).Unit = 'rad';
model_train_initial.InitialStates(2).Name = 'theta_2'; model_train_initial.InitialStates(2).Unit = 'rad';
model_train_initial.InitialStates(3).Name = 'theta_1_dot'; model_train_initial.InitialStates(3).Unit = 'rad/s';
model_train_initial.InitialStates(4).Name = 'theta_2_dot'; model_train_initial.InitialStates(4).Unit = 'rad/s';

% Fix initial conditions for training
model_train_initial = setinit(model_train_initial, 'Value', initial_states_train);
model_train_initial = setinit(model_train_initial, 'Fixed', {true; true; true; true});

% --- 5. Estimate Parameters ---
opt = nlgreyestOptions;
opt.Display = 'Full';
opt.EstimateCovariance = true;
opt.SearchOptions.MaxIterations = 150;

% [UPDATE]: Ensure Levenberg-Marquardt is set correctly on the main object
opt.SearchMethod = 'lm';

% OutputWeight prioritizes the passive link slightly to ensure it doesn't destabilize
opt.OutputWeight = diag([1, 8]);

fprintf('\nStarting grey-box estimation for active linkage (Levenberg-Marquardt)...\n');
model_est = nlgreyest(iddata_train, model_train_initial, opt);

% --- 6. Validate & Plot ---
model_val = model_est;
model_val = setinit(model_val, 'Value', initial_states_val);
model_val = setinit(model_val, 'Fixed', {true; true; true; true});

compare_opt = compareOptions;
compare_opt.InitialCondition = 'estimate';

[y_train_compare, fit_train] = compare(iddata_train, model_est, compare_opt);
[y_val_compare, fit_val]     = compare(iddata_val, model_val, compare_opt);

fprintf('\nTraining fit (theta1, theta2):\n'); disp(fit_train);
fprintf('\nValidation fit (theta1, theta2):\n'); disp(fit_val);

figure('Name', 'Training comparison', 'Units', 'normalized', 'Position', [0.05, 0.08, 0.85, 0.8]);
compare(iddata_train, model_est, compare_opt);
title('Training comparison');

figure('Name', 'Validation comparison (t=60s to end)', 'Units', 'normalized', 'Position', [0.08, 0.08, 0.85, 0.8]);
compare(iddata_val, model_val, compare_opt);
title('Validation comparison (Truncated)');

%% --- Print Detailed Estimation Results ---
fprintf('\n======================================================\n');
fprintf('             FINAL ESTIMATION STATISTICS                \n');
fprintf('======================================================\n');

num_params = length(model_est.Parameters);
fprintf('\nEstimated Active Parameters (± 1 Standard Deviation):\n');
fprintf('------------------------------------------------------\n');

try
    cov_matrix = getcov(model_est);
    std_devs = sqrt(diag(cov_matrix));
    cov_idx = 1;
    
    for k = 1:num_params
        if ~model_est.Parameters(k).Fixed
            fprintf('%-8s : %10.4f  ±  %-10.4f\n', ...
                model_est.Parameters(k).Name, ...
                model_est.Parameters(k).Value, ...
                std_devs(cov_idx));
            cov_idx = cov_idx + 1;
        end
    end
catch
    fprintf('Covariance could not be estimated. Showing values only:\n');
    for k = 1:6
        fprintf('%-8s : %10.4f\n', model_est.Parameters(k).Name, model_est.Parameters(k).Value);
    end
end

fprintf('\nFinal Loss Function (Cost): %g\n', model_est.Report.Fit.LossFcn);
fprintf('======================================================\n\n');

% --- 7. Save Models and Stop Logging ---
save('active_rod_identification_result.mat', 'model_est', 'model_val');
fprintf('Model saved to active_rod_identification_result.mat\n');
diary off;