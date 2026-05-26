clear;
clear functions;
clc;
close all;

% Load measurement data
data_train = load("identification_data\prbs_train_121s_amp_015.mat");
data_val = load("identification_data\prbs_validation_121s_amp_015.mat");

% Define sampling period
Ts = 0.01;

% Extract data
t_train = data_train.theta_1.Time(:);
theta_1_train = data_train.theta_1.Data(:);
theta_2_train = data_train.theta_2.Data(:);
u_train = data_train.u_ts.Data(1:12001);

t_val = data_val.theta_1.Time(:);
theta_1_val = data_val.theta_1.Data(:);
theta_2_val = data_val.theta_2.Data(:);
u_val = data_val.u_ts.Data(1:12001);

%% Make sure all signals have equal length
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

% Plot training data
figure('Name', 'Training data', 'Units', 'normalized', 'Position', [0.05, 0.08, 0.85, 0.8]);

tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t_train, theta_1_train, 'LineWidth', 1.1);
grid on;
ylabel('\theta_1 [rad]');
title('Training data');

nexttile;
plot(t_train, theta_2_train, 'LineWidth', 1.1);
grid on;
ylabel('\theta_2 [rad]');

nexttile;
plot(t_train, u_train, 'LineWidth', 1.1);
grid on;
ylabel('u');
xlabel('Time [s]');

% Plot validation data
figure('Name', 'Validation data', 'Units', 'normalized', 'Position', [0.08, 0.08, 0.85, 0.8]);

tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t_val, theta_1_val, 'LineWidth', 1.1);
grid on;
ylabel('\theta_1 [rad]');
title('Validation data');

nexttile;
plot(t_val, theta_2_val, 'LineWidth', 1.1);
grid on;
ylabel('\theta_2 [rad]');

nexttile;
plot(t_val, u_val, 'LineWidth', 1.1);
grid on;
ylabel('u');
xlabel('Time [s]');

% Plot fourier spectra
fs = 1 / Ts;              % Sampling frequency [Hz]
maxFreqToShow = 10;       % Show spectrum up to 20 Hz


% Training spectrum: theta_1
signal_train_theta1 = theta_1_train(:);
signal_train_theta1 = signal_train_theta1 - mean(signal_train_theta1, 'omitnan');

N_train_theta1 = length(signal_train_theta1);
Y_train_theta1 = fft(signal_train_theta1);

P2_train_theta1 = abs(Y_train_theta1 / N_train_theta1);
P_train_theta1 = P2_train_theta1(1:floor(N_train_theta1/2) + 1);

if length(P_train_theta1) > 2
    P_train_theta1(2:end-1) = 2 * P_train_theta1(2:end-1);
end

f_train_theta1 = fs * (0:floor(N_train_theta1/2))' / N_train_theta1;


% Training spectrum: theta_2
signal_train_theta2 = theta_2_train(:);
signal_train_theta2 = signal_train_theta2 - mean(signal_train_theta2, 'omitnan');

N_train_theta2 = length(signal_train_theta2);
Y_train_theta2 = fft(signal_train_theta2);

P2_train_theta2 = abs(Y_train_theta2 / N_train_theta2);
P_train_theta2 = P2_train_theta2(1:floor(N_train_theta2/2) + 1);

if length(P_train_theta2) > 2
    P_train_theta2(2:end-1) = 2 * P_train_theta2(2:end-1);
end

f_train_theta2 = fs * (0:floor(N_train_theta2/2))' / N_train_theta2;


% Training spectrum: input u
signal_train_u = u_train(:);
signal_train_u = signal_train_u - mean(signal_train_u, 'omitnan');

N_train_u = length(signal_train_u);
Y_train_u = fft(signal_train_u);

P2_train_u = abs(Y_train_u / N_train_u);
P_train_u = P2_train_u(1:floor(N_train_u/2) + 1);

if length(P_train_u) > 2
    P_train_u(2:end-1) = 2 * P_train_u(2:end-1);
end

f_train_u = fs * (0:floor(N_train_u/2))' / N_train_u;


% Plot Fourier spectra of training data
figure('Name', 'Training data Fourier spectra', ...
       'Units', 'normalized', ...
       'Position', [0.05, 0.08, 0.85, 0.8]);

tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(f_train_theta1, P_train_theta1, 'LineWidth', 1.1);
grid on;
ylabel('Amplitude');
title('Training Fourier spectra');
xlim([0, maxFreqToShow]);
legend('\theta_1', 'Location', 'best');

nexttile;
plot(f_train_theta2, P_train_theta2, 'LineWidth', 1.1);
grid on;
ylabel('Amplitude');
xlim([0, maxFreqToShow]);
legend('\theta_2', 'Location', 'best');

nexttile;
plot(f_train_u, P_train_u, 'LineWidth', 1.1);
grid on;
ylabel('Amplitude');
xlabel('Frequency [Hz]');
xlim([0, maxFreqToShow]);
legend('u', 'Location', 'best');


% Validation spectrum: theta_1
signal_val_theta1 = theta_1_val(:);
signal_val_theta1 = signal_val_theta1 - mean(signal_val_theta1, 'omitnan');

N_val_theta1 = length(signal_val_theta1);
Y_val_theta1 = fft(signal_val_theta1);

P2_val_theta1 = abs(Y_val_theta1 / N_val_theta1);
P_val_theta1 = P2_val_theta1(1:floor(N_val_theta1/2) + 1);

if length(P_val_theta1) > 2
    P_val_theta1(2:end-1) = 2 * P_val_theta1(2:end-1);
end

f_val_theta1 = fs * (0:floor(N_val_theta1/2))' / N_val_theta1;


% Validation spectrum: theta_2
signal_val_theta2 = theta_2_val(:);
signal_val_theta2 = signal_val_theta2 - mean(signal_val_theta2, 'omitnan');

N_val_theta2 = length(signal_val_theta2);
Y_val_theta2 = fft(signal_val_theta2);

P2_val_theta2 = abs(Y_val_theta2 / N_val_theta2);
P_val_theta2 = P2_val_theta2(1:floor(N_val_theta2/2) + 1);

if length(P_val_theta2) > 2
    P_val_theta2(2:end-1) = 2 * P_val_theta2(2:end-1);
end

f_val_theta2 = fs * (0:floor(N_val_theta2/2))' / N_val_theta2;


% Validation spectrum: input u
signal_val_u = u_val(:);
signal_val_u = signal_val_u - mean(signal_val_u, 'omitnan');

N_val_u = length(signal_val_u);
Y_val_u = fft(signal_val_u);

P2_val_u = abs(Y_val_u / N_val_u);
P_val_u = P2_val_u(1:floor(N_val_u/2) + 1);

if length(P_val_u) > 2
    P_val_u(2:end-1) = 2 * P_val_u(2:end-1);
end

f_val_u = fs * (0:floor(N_val_u/2))' / N_val_u;


% Plot Fourier spectra of validation data
figure('Name', 'Validation data Fourier spectra', ...
       'Units', 'normalized', ...
       'Position', [0.08, 0.08, 0.85, 0.8]);

tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(f_val_theta1, P_val_theta1, 'LineWidth', 1.1);
grid on;
ylabel('Amplitude');
title('Validation Fourier spectra');
xlim([0, maxFreqToShow]);
legend('\theta_1', 'Location', 'best');

nexttile;
plot(f_val_theta2, P_val_theta2, 'LineWidth', 1.1);
grid on;
ylabel('Amplitude');
xlim([0, maxFreqToShow]);
legend('\theta_2', 'Location', 'best');

nexttile;
plot(f_val_u, P_val_u, 'LineWidth', 1.1);
grid on;
ylabel('Amplitude');
xlabel('Frequency [Hz]');
xlim([0, maxFreqToShow]);
legend('u', 'Location', 'best');

% Load data into iddata
y_train = [theta_1_train, theta_2_train];
y_val = [theta_1_val, theta_2_val];

iddata_train = iddata(y_train, u_train, Ts);
iddata_val = iddata(y_val, u_val, Ts);

iddata_train.InputName  = {'u'};
iddata_train.OutputName = {'theta_1', 'theta_2'};
iddata_train.InputUnit  = {'V'};
iddata_train.OutputUnit = {'rad', 'rad'};
iddata_train.TimeUnit   = 's';

iddata_val.InputName  = {'u'};
iddata_val.OutputName = {'theta_1', 'theta_2'};
iddata_val.InputUnit  = {'V'};
iddata_val.OutputUnit = {'rad', 'rad'};
iddata_val.TimeUnit   = 's';

% Set initial conditions

% Since your model uses relative coordinates:
% x(1) = theta_1
% x(2) = theta_2
% x(3) = theta_1_dot
% x(4) = theta_2_dot

initial_states_train = {
    theta_1_train(1)
    theta_2_train(1)
    0
    0
};

initial_states_val = {
    theta_1_val(1)
    theta_2_val(1)
    0
    0
};


% Create nonlinear grey-box model
model_file = 'greybox_id_6par_model';
order = [2 1 4];     % 2 outputs, 1 input, 4 states
Ts_model = 0;        % continuous-time model


% Set initial parameter guesses
p_a_initial  = 14;
p_b1_initial = 670;
p_b2_initial = 0.5;
p_g1_initial = 97;
p_g2_initial = 102;
p_u_initial  = 4000;

parameters = {
    p_a_initial
    p_b1_initial
    p_b2_initial
    p_g1_initial
    p_g2_initial
    p_u_initial
};


% Create initial model
model_train_initial = idnlgrey(model_file, order, parameters, initial_states_train, Ts_model);


% Configure parameter names and bounds
parameter_names = {
    'p_a'
    'p_b1'
    'p_b2'
    'p_g1'
    'p_g2'
    'p_u'
};

% Minimum bounds (Enforcing physical reality)
minimum_values = [
    5       % p_a  (Strictly positive to prevent matrix singularity)
    100     % p_b1 (Motor back-EMF guarantees substantial damping)
    0       % p_b2 (Passive bearing friction cannot be negative)
    20      % p_g1 (Rod 1 has measurable mass)
    20      % p_g2 (Rod 2 has measurable mass)
    500     % p_u  (Assuming correct motor wiring, gain must be positive)
];

% Maximum bounds (Allowing reasonable search space)
maximum_values = [
    40      % p_a  (Allows for heavier rods/hubs than estimated)
    1500    % p_b1 (Allows for a highly resistive motor)
    20      % p_b2 (Passive joint damping should remain very small)
    200     % p_g1 (Allows for up to ~2x the estimated mass)
    200     % p_g2 (Allows for up to ~2x the estimated mass)
    8000    % p_u  (Allows for a highly powerful motor/amplifier combo)
];

for k = 1:length(parameter_names)
    model_train_initial.Parameters(k).Name = parameter_names{k};
    model_train_initial.Parameters(k).Minimum = minimum_values(k);
    model_train_initial.Parameters(k).Maximum = maximum_values(k);
    model_train_initial.Parameters(k).Fixed = false;
end


% Configure input/output/state names

model_train_initial.InputName = {'u'};
model_train_initial.InputUnit = {'V'};

model_train_initial.OutputName = {'theta_1', 'theta_2'};
model_train_initial.OutputUnit = {'rad', 'rad'};

model_train_initial.TimeUnit = 's';

model_train_initial.InitialStates(1).Name = 'theta_1';
model_train_initial.InitialStates(1).Unit = 'rad';

model_train_initial.InitialStates(2).Name = 'theta_2';
model_train_initial.InitialStates(2).Unit = 'rad';

model_train_initial.InitialStates(3).Name = 'theta_1_dot';
model_train_initial.InitialStates(3).Unit = 'rad/s';

model_train_initial.InitialStates(4).Name = 'theta_2_dot';
model_train_initial.InitialStates(4).Unit = 'rad/s';


% Fix initial conditions
model_train_initial = setinit(model_train_initial, 'Value', initial_states_train);
model_train_initial = setinit(model_train_initial, 'Fixed', {true; true; true; true});


%Set estimator options
opt = nlgreyestOptions;
opt.Display = 'Full';
opt.EstimateCovariance = true;
opt.SearchOptions.MaxIterations = 150;

% Optional weighting. Use this if theta_2 is otherwise fitted poorly.
opt.OutputWeight = diag([1, 8]);


% Estimate model on training data
model_est = nlgreyest(iddata_train, model_train_initial, opt);


% Create validation model using same estimated parameters but validation initial state
model_val = model_est;

model_val = setinit(model_val, 'Value', initial_states_val);
model_val = setinit(model_val, 'Fixed', {true; true; true; true});


% Compare training and validation
compare_opt = compareOptions;
compare_opt.InitialCondition = 'model';

[y_train_compare, fit_train] = compare(iddata_train, model_est, compare_opt);
[y_val_compare, fit_val]     = compare(iddata_val, model_val, compare_opt);

fprintf('\nTraining fit:\n');
disp(fit_train);

fprintf('\nValidation fit:\n');
disp(fit_val);


% Plot comparisons
figure('Name', 'Training comparison', 'Units', 'normalized', 'Position', [0.05, 0.08, 0.85, 0.8]);

compare(iddata_train, model_est, compare_opt);
title('Training comparison');


figure('Name', 'Validation comparison', 'Units', 'normalized', 'Position', [0.08, 0.08, 0.85, 0.8]);

compare(iddata_val, model_val, compare_opt);
title('Validation comparison');

%% --- Print Detailed Estimation Results ---
fprintf('\n======================================================\n');
fprintf('             FINAL ESTIMATION STATISTICS                \n');
fprintf('======================================================\n');

% 1. Extract Parameter Names and Final Values
num_params = length(model_est.Parameters);
param_names = cell(num_params, 1);
param_vals = zeros(num_params, 1);

for k = 1:num_params
    param_names{k} = model_est.Parameters(k).Name;
    param_vals(k) = model_est.Parameters(k).Value;
end

% 2. Extract Covariance and Calculate Standard Deviations
% getcov() returns the covariance matrix for all FREE parameters.
% Since you set Fixed = false for all 6 parameters, this will be a 6x6 matrix.
cov_matrix = getcov(model_est);
std_devs = sqrt(diag(cov_matrix));

% 3. Print the Parameter Values with their Uncertainties (± 1 Std Dev)
fprintf('\nEstimated Parameters (± 1 Standard Deviation):\n');
fprintf('------------------------------------------------------\n');
for k = 1:num_params
    % Formatted to align nicely in the console
    fprintf('%-8s : %10.4f  ±  %-10.4f\n', param_names{k}, param_vals(k), std_devs(k));
end

% 4. Print the Full Covariance Matrix
fprintf('\nCovariance Matrix (6x6):\n');
fprintf('------------------------------------------------------\n');
disp(cov_matrix);

% 5. Print the Loss Function (Cost)
% This represents the final value of the determinant of the prediction error covariance
fprintf('\nFinal Loss Function (Cost): %g\n', model_est.Report.Fit.LossFcn);
fprintf('======================================================\n\n');
