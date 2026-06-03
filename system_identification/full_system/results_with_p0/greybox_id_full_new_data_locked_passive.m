clear;
clear functions;
clc;
close all;

%% ================================================================
% Full-system nonlinear grey-box identification using new data
%
% This script uses the cleaned 30 s PRBS and chirp measurements.
%
% Expected folders:
%   prbs/
%       fullsystem_prbs_A0p16_run01.mat
%       ...
%       fullsystem_prbs_A0p34_run10.mat
%
%   chirp/
%       fullsystem_chirp_A0p16_run01.mat
%       ...
%       fullsystem_chirp_A0p34_run10.mat
%
% Each file must contain:
%   theta_1
%   theta_2
%   u_ts
%
% Identification split:
%   odd runs  = identification data
%   even runs = validation data
%
% Passive-link parameters are fixed to the values found during the
% passive-link identification experiment. Only the active/full-system
% parameters are estimated.
% ================================================================

diary('full_system_id_new_data_console_output.txt');

fprintf('\n======================================================\n');
fprintf('Starting full-system identification with new data\n');
fprintf('Passive-link parameters are locked/fixed\n');
fprintf('======================================================\n\n');

%% Settings

Ts = 0.01;
inputSign = -1;     % old script used u = -u. Try +1 if the fit is clearly mirrored.

amplitudes = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];

idxEst = [1 3 5 7 9];      % odd runs for identification
idxVal = [2 4 6 8 10];     % even runs for validation

fprintf('Identification runs: 1, 3, 5, 7, 9 for both PRBS and chirp\n');
fprintf('Validation runs:     2, 4, 6, 8, 10 for both PRBS and chirp\n\n');

%% ================================================================
% Load data and build multi-experiment iddata objects
% ================================================================

zEst = [];
zVal = [];

x0Est = [];
x0Val = [];

nEst = 0;
nVal = 0;

for k = 1:10

    ampText = strrep(sprintf('%.2f', amplitudes(k)), '.', 'p');

    prbsFile  = fullfile('measurement_data\prbs',  sprintf('fullsystem_prbs_A%s_run%02d.mat',  ampText, k));
    chirpFile = fullfile('measurement_data\chirp', sprintf('fullsystem_chirp_A%s_run%02d.mat', ampText, k));

    %% PRBS file
    load(prbsFile, 'theta_1', 'theta_2', 'u_ts');

    theta1 = double(squeeze(theta_1.Data(:)));
    theta2 = double(squeeze(theta_2.Data(:)));
    u = inputSign * double(squeeze(u_ts.Data(:)));

    y = [theta1, theta2];
    z = iddata(y, u, Ts);
    z.Name = sprintf('prbs_run_%02d', k);
    z.InputName = {'u'};
    z.InputUnit = {'V'};
    z.OutputName = {'theta_1', 'theta_2'};
    z.OutputUnit = {'rad', 'rad'};
    z.TimeUnit = 's';

    x0 = [theta1(1); theta2(1); 0; 0];

    if ismember(k, idxEst)
        nEst = nEst + 1;
        if nEst == 1
            zEst = z;
        else
            zEst = merge(zEst, z);
        end
        x0Est(:, nEst) = x0;
    else
        nVal = nVal + 1;
        if nVal == 1
            zVal = z;
        else
            zVal = merge(zVal, z);
        end
        x0Val(:, nVal) = x0;
    end

    %% Chirp file
    load(chirpFile, 'theta_1', 'theta_2', 'u_ts');

    theta1 = double(squeeze(theta_1.Data(:)));
    theta2 = double(squeeze(theta_2.Data(:)));
    u = inputSign * double(squeeze(u_ts.Data(:)));

    y = [theta1, theta2];
    z = iddata(y, u, Ts);
    z.Name = sprintf('chirp_run_%02d', k);
    z.InputName = {'u'};
    z.InputUnit = {'V'};
    z.OutputName = {'theta_1', 'theta_2'};
    z.OutputUnit = {'rad', 'rad'};
    z.TimeUnit = 's';

    x0 = [theta1(1); theta2(1); 0; 0];

    if ismember(k, idxEst)
        nEst = nEst + 1;
        zEst = merge(zEst, z);
        x0Est(:, nEst) = x0;
    else
        nVal = nVal + 1;
        zVal = merge(zVal, z);
        x0Val(:, nVal) = x0;
    end
end

fprintf('Number of identification experiments: %d\n', nEst);
fprintf('Number of validation experiments:     %d\n\n', nVal);

%% ================================================================
% Passive-link parameters from passive-link identification
% These are locked during full-system identification.
% ================================================================

p_b2_val      = 0.067028846;
p_g2_val      = 112.02675;
p_c2_val      = 0.24908901;
p_sdelta2_val = 0.18360238;
v_s2_val      = 3.5;
eps_v2_val    = 0.03;

%% ================================================================
% Grey-box model setup
% ================================================================

model_file = 'greybox_id_full_model';
order = [2 1 4];       % 2 outputs, 1 input, 4 states
Ts_model = 0;          % continuous-time grey-box model

% Parameters:
%  1 p_a        active/full-system inertia ratio
%  2 p_b1       active joint effective viscous damping
%  3 p_c1       active joint Coulomb friction
%  4 p_g1       active-link gravity parameter
%  5 p_u        motor input gain
%  6 p_0        constant torque/input bias
%  7 p_b2       passive joint viscous damping, fixed
%  8 p_g2       passive-link gravity parameter, fixed
%  9 p_c2       passive joint Coulomb friction, fixed
% 10 p_sdelta2  passive Stribeck static-minus-Coulomb term, fixed
% 11 v_s2       passive Stribeck velocity, fixed
% 12 eps_v1     active smoothing velocity, fixed
% 13 eps_v2     passive smoothing velocity, fixed

parameters = {
    14;              % 1: p_a
    670;             % 2: p_b1
    0;               % 3: p_c1, fixed to zero in stage 1
    97;              % 4: p_g1
    4000;            % 5: p_u
    0;               % 6: p_0
    p_b2_val;        % 7: p_b2, fixed
    p_g2_val;        % 8: p_g2, fixed
    p_c2_val;        % 9: p_c2, fixed
    p_sdelta2_val;   % 10: p_sdelta2, fixed
    v_s2_val;        % 11: v_s2, fixed
    0.05;            % 12: eps_v1, fixed
    eps_v2_val       % 13: eps_v2, fixed
};

initial_states_est = {
    x0Est(1, :);
    x0Est(2, :);
    x0Est(3, :);
    x0Est(4, :)
};

model0 = idnlgrey(model_file, order, parameters, initial_states_est, Ts_model);

parameter_names = {'p_a', 'p_b1', 'p_c1', 'p_g1', 'p_u', 'p_0', ...
    'p_b2', 'p_g2', 'p_c2', 'p_sdelta2', 'v_s2', 'eps_v1', 'eps_v2'};

minimum_values = [5, 100, 0, 5, 500, -200, 0, 0, 0, 0, 0.02, 0.001, 0.001];
maximum_values = [80, 4000, 400, 400, 30000, 200, Inf, Inf, Inf, Inf, 20, 1, 1];

% Stage 1:
% Estimate p_a, p_b1, p_g1, p_u and p_0.
% Keep p_c1 fixed at zero first, and keep all passive-link parameters fixed.
fixed_stage1 = [false, false, true, false, false, false, ...
                true,  true,  true, true,  true,  true, true];

for k = 1:length(parameter_names)
    model0.Parameters(k).Name = parameter_names{k};
    model0.Parameters(k).Minimum = minimum_values(k);
    model0.Parameters(k).Maximum = maximum_values(k);
    model0.Parameters(k).Fixed = fixed_stage1(k);
end

model0.InputName = {'u'};
model0.InputUnit = {'V'};
model0.OutputName = {'theta_1', 'theta_2'};
model0.OutputUnit = {'rad', 'rad'};
model0.TimeUnit = 's';

model0.InitialStates(1).Name = 'theta_1';
model0.InitialStates(1).Unit = 'rad';
model0.InitialStates(2).Name = 'theta_2';
model0.InitialStates(2).Unit = 'rad';
model0.InitialStates(3).Name = 'theta_1_dot';
model0.InitialStates(3).Unit = 'rad/s';
model0.InitialStates(4).Name = 'theta_2_dot';
model0.InitialStates(4).Unit = 'rad/s';

model0 = setinit(model0, 'Fixed', {
    true(1, nEst);
    true(1, nEst);
    true(1, nEst);
    true(1, nEst)
});

%% Compare initial model before estimation

compare_opt = compareOptions;
compare_opt.InitialCondition = 'estimate';

figure('Name', 'Initial model fit before estimation', 'Units', 'normalized', 'Position', [0.05 0.08 0.88 0.78]);
compare(zEst, model0, compare_opt);
title('Initial model fit before estimation');

%% ================================================================
% Estimation options
% ================================================================

opt1 = nlgreyestOptions;
opt1.Display = 'Full';
opt1.EstimateCovariance = true;
opt1.SearchMethod = 'lm';
opt1.SearchOptions.MaxIterations = 100;
opt1.OutputWeight = diag([1, 8]);

opt2 = nlgreyestOptions;
opt2.Display = 'Full';
opt2.EstimateCovariance = true;
opt2.SearchMethod = 'lm';
opt2.SearchOptions.MaxIterations = 150;
opt2.OutputWeight = diag([1, 8]);

%% ================================================================
% Stage 1: estimate active linear-ish parameters with p_c1 fixed
% ================================================================

fprintf('\n======================================================\n');
fprintf('Stage 1: estimating active parameters with p_c1 fixed\n');
fprintf('Free parameters: p_a, p_b1, p_g1, p_u, p_0\n');
fprintf('Fixed parameters: p_c1 and all passive-link parameters\n');
fprintf('======================================================\n\n');

model_stage1 = nlgreyest(zEst, model0, opt1);
model_stage1.Name = 'Stage 1 active model, p_c1 fixed';

fprintf('\nStage 1 parameter values:\n');
for k = 1:length(model_stage1.Parameters)
    if model_stage1.Parameters(k).Fixed
        fprintf('%-10s = %12.6g   fixed\n', model_stage1.Parameters(k).Name, model_stage1.Parameters(k).Value);
    else
        fprintf('%-10s = %12.6g   estimated\n', model_stage1.Parameters(k).Name, model_stage1.Parameters(k).Value);
    end
end

[~, fitStage1Est] = compare(zEst, model_stage1, compare_opt);

fitStage1Values = [];
if iscell(fitStage1Est)
    for k = 1:numel(fitStage1Est)
        fitStage1Values = [fitStage1Values; fitStage1Est{k}(:)]; %#ok<AGROW>
    end
else
    fitStage1Values = fitStage1Est(:);
end
fitStage1Values = fitStage1Values(isfinite(fitStage1Values));

fprintf('\nStage 1 mean identification fit: %.2f %%\n', mean(fitStage1Values));

figure('Name', 'Stage 1 identification fit', 'Units', 'normalized', 'Position', [0.05 0.08 0.88 0.78]);
compare(zEst, model_stage1, compare_opt);
title('Stage 1 identification fit, p_{c1} fixed');

%% ================================================================
% Stage 2: unlock p_c1 and estimate final active parameter set
% ================================================================

fprintf('\n======================================================\n');
fprintf('Stage 2: estimating final active parameter set\n');
fprintf('Free parameters: p_a, p_b1, p_c1, p_g1, p_u, p_0\n');
fprintf('Fixed parameters: all passive-link parameters and smoothing constants\n');
fprintf('======================================================\n\n');

model_stage2_start = model_stage1;
model_stage2_start.Parameters(3).Value = 10;
model_stage2_start.Parameters(3).Fixed = false;

% Make extra sure the passive-link parameters remain locked.
for k = 7:13
    model_stage2_start.Parameters(k).Fixed = true;
end

model_est = nlgreyest(zEst, model_stage2_start, opt2);
model_est.Name = 'Final full-system model with locked passive-link parameters';

%% ================================================================
% Validation model: same parameters, validation initial states
% ================================================================

initial_states_val = {
    x0Val(1, :);
    x0Val(2, :);
    x0Val(3, :);
    x0Val(4, :)
};

model_val = model_est;
model_val = setinit(model_val, 'Value', initial_states_val);
model_val = setinit(model_val, 'Fixed', {
    true(1, nVal);
    true(1, nVal);
    true(1, nVal);
    true(1, nVal)
});

%% ================================================================
% Compare final model on identification and validation data
% ================================================================

[yEst, fitEst] = compare(zEst, model_est, compare_opt);
[yVal, fitVal] = compare(zVal, model_val, compare_opt);

fitEstValues = [];
if iscell(fitEst)
    for k = 1:numel(fitEst)
        fitEstValues = [fitEstValues; fitEst{k}(:)]; %#ok<AGROW>
    end
else
    fitEstValues = fitEst(:);
end
fitEstValues = fitEstValues(isfinite(fitEstValues));

fitValValues = [];
if iscell(fitVal)
    for k = 1:numel(fitVal)
        fitValValues = [fitValValues; fitVal{k}(:)]; %#ok<AGROW>
    end
else
    fitValValues = fitVal(:);
end
fitValValues = fitValValues(isfinite(fitValValues));

fprintf('\nFinal identification fit, per experiment and output:\n');
disp(fitEst);

fprintf('\nFinal validation fit, per experiment and output:\n');
disp(fitVal);

fprintf('\nFinal mean identification fit: %.2f %%\n', mean(fitEstValues));
fprintf('Final mean validation fit:     %.2f %%\n', mean(fitValValues));

figure('Name', 'Final identification fit', 'Units', 'normalized', 'Position', [0.05 0.08 0.88 0.78]);
compare(zEst, model_est, compare_opt);
title('Final identification fit');

figure('Name', 'Final validation fit', 'Units', 'normalized', 'Position', [0.08 0.08 0.88 0.78]);
compare(zVal, model_val, compare_opt);
title('Final validation fit');

%% ================================================================
% Print parameter values, covariance, standard deviations and correlations
% ================================================================

fprintf('\n======================================================\n');
fprintf('Final parameter values and uncertainty\n');
fprintf('======================================================\n\n');

nPar = length(model_est.Parameters);
parameter = cell(nPar, 1);
value = zeros(nPar, 1);
fixed = false(nPar, 1);
variance = NaN(nPar, 1);
stdDev = NaN(nPar, 1);
relativeStdPercent = NaN(nPar, 1);

for k = 1:nPar
    parameter{k} = model_est.Parameters(k).Name;
    value(k) = model_est.Parameters(k).Value;
    fixed(k) = model_est.Parameters(k).Fixed;
end

freeIdx = find(~fixed);
covFree = getcov(model_est, 'value', 'free');

for k = 1:length(freeIdx)
    variance(freeIdx(k)) = covFree(k, k);
    stdDev(freeIdx(k)) = sqrt(max(covFree(k, k), 0));
    relativeStdPercent(freeIdx(k)) = 100 * stdDev(freeIdx(k)) / abs(value(freeIdx(k)));
end

resultTable = table(parameter, value, fixed, variance, stdDev, relativeStdPercent, ...
    'VariableNames', {'Parameter', 'Value', 'Fixed', 'Variance', 'StdDev', 'RelativeStdPercent'});

disp(resultTable);

fprintf('\nFree-parameter covariance matrix:\n');
disp(covFree);

freeNames = parameter(freeIdx);
stdFree = sqrt(max(diag(covFree), 0));
corrFree = covFree ./ (stdFree * stdFree.');
corrFree(1:size(corrFree, 1)+1:end) = 1;

corrVarNames = matlab.lang.makeValidName(freeNames);
corrTable = array2table(corrFree, 'VariableNames', corrVarNames, 'RowNames', freeNames);

fprintf('\nFree-parameter correlation matrix:\n');
disp(corrTable);

fprintf('\nFinal loss function: %g\n', model_est.Report.Fit.LossFcn);

%% ================================================================
% Save results
% ================================================================

save('full_system_id_new_data_locked_passive_result.mat', ...
    'model_est', 'model_val', 'model_stage1', 'zEst', 'zVal', ...
    'fitEst', 'fitVal', 'fitStage1Est', 'resultTable', 'covFree', 'corrTable', ...
    'idxEst', 'idxVal', 'amplitudes', 'inputSign');

writetable(resultTable, 'full_system_id_new_data_locked_passive_parameters.csv');
writetable(corrTable, 'full_system_id_new_data_locked_passive_correlations.csv', 'WriteRowNames', true);

fprintf('\nSaved result file: full_system_id_new_data_locked_passive_result.mat\n');
fprintf('Saved parameter table: full_system_id_new_data_locked_passive_parameters.csv\n');
fprintf('Saved correlation table: full_system_id_new_data_locked_passive_correlations.csv\n');

diary off;
