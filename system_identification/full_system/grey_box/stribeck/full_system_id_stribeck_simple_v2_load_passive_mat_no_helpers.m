clear;
clear functions;
clc;
close all;

%% ================================================================
% Full-system nonlinear grey-box identification, Stribeck friction version
%
% This is a cleaned version of the original full-system identification
% script. The identification functionality is kept the same:
%   - cleaned PRBS and chirp data are used
%   - odd runs are used for identification
%   - even runs are used for validation
%   - passive-link parameters are fixed
%   - active/full-system parameters are estimated in two stages
%   - stage 1 keeps p_c1 fixed
%   - stage 2 unlocks p_c1
%   - compare uses estimated initial conditions
%
% Required model file on the MATLAB path:
%   greybox_id_full_stribeck_model.m
% ================================================================

%% ================================================================
% Configuration
% ================================================================

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end

% Change these paths if your folder structure is different.
dataFolder   = fullfile(scriptFolder, 'measurement_data');
prbsFolder   = fullfile(dataFolder, 'prbs');
chirpFolder  = fullfile(dataFolder, 'chirp');
modelsFolder = scriptFolder;
outputFolder = scriptFolder;

% Example if your files are stored one folder above the script:
% dataFolder   = fullfile(scriptFolder, '..', '..', 'measurement_data');
% prbsFolder   = fullfile(dataFolder, 'prbs');
% chirpFolder  = fullfile(dataFolder, 'chirp');
% modelsFolder = fullfile(scriptFolder, '..', 'models');

addpath(modelsFolder);

modelFile = 'greybox_id_full_stribeck_model';

Ts = 0.01;
TsModel = 0;         % continuous-time grey-box model
inputSign = -1;      % original script used u = -u. Try +1 if the fit is mirrored.

amplitudes = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];

idxEst = [1 3 5 7 9];    % odd runs for identification
idxVal = [2 4 6 8 10];   % even runs for validation

runPrbs = true;
runChirp = true;

useDiary = true;
consoleLogFile = 'full_system_id_stribeck_console_output.txt';

showInitialCompare = true;
showStage1Compare = true;
showFinalCompare = true;

figureUnits = 'pixels';
figurePosition = [80 80 1150 720];

% Passive-link parameters are loaded from the passive-link ID result file.
% Keep the manual values as a fallback only.
loadPassiveParametersFromMat = true;
passiveResultMatFile = fullfile(scriptFolder, 'passive_link_3step_stribeck_result_v5_no_helpers.mat');

% Manual fallback values. These are only used when
% loadPassiveParametersFromMat = false.
p_b2_val_manual      = 0.067028846;
p_g2_val_manual      = 112.02675;
p_c2_val_manual      = 0.24908901;
p_sdelta2_val_manual = 0.18360238;
v_s2_val_manual      = 3.5;
eps_v2_val_manual    = 0.03;

% These variables are assigned after the model-file check below.
p_b2_val      = p_b2_val_manual;
p_g2_val      = p_g2_val_manual;
p_c2_val      = p_c2_val_manual;
p_sdelta2_val = p_sdelta2_val_manual;
v_s2_val      = v_s2_val_manual;
eps_v2_val    = eps_v2_val_manual;

% Initial full-system parameter guesses.
p_a_init   = 14;
p_b1_init  = 670;
p_c1_init  = 0;
p_g1_init  = 97;
p_u_init   = 4000;
p_0_init   = 0;
eps_v1_val = 0.05;

% Bounds for the Stribeck full-system model parameters:
% [p_a, p_b1, p_c1, p_g1, p_u, p_0, p_b2, p_g2, p_c2, p_sdelta2, v_s2, eps_v1, eps_v2]
minimumValues = [5, 100, 0, 5, 500, -200, 0, 0, 0, 0, 0.02, 0.001, 0.001];
maximumValues = [80, 4000, 800, 500, 30000, 200, Inf, Inf, Inf, Inf, 20, 1, 1];

% Stage 1:
% Estimate p_a, p_b1, p_g1, p_u and p_0.
% Keep p_c1 fixed at zero first, and keep all passive-link parameters fixed.
fixedStage1 = [false, false, true, false, false, false, ...
               true,  true,  true, true,  true,  true, true];

% Stage 2:
% Unlock p_c1 and estimate final active/full-system parameter set.
fixedStage2 = [false, false, false, false, false, false, ...
               true,  true,  true,  true,  true,  true, true];

stage2_p_c1_initial_value = 10;

% Estimation options.
optStage1 = nlgreyestOptions;
optStage1.Display = 'Full';
optStage1.EstimateCovariance = true;
optStage1.SearchMethod = 'lm';
optStage1.SearchOptions.MaxIterations = 100;
optStage1.OutputWeight = diag([1, 8]);

optStage2 = nlgreyestOptions;
optStage2.Display = 'Full';
optStage2.EstimateCovariance = true;
optStage2.SearchMethod = 'lm';
optStage2.SearchOptions.MaxIterations = 150;
optStage2.OutputWeight = diag([1, 8]);

compareOpt = compareOptions;
compareOpt.InitialCondition = 'estimate';

resultMatFile = fullfile(outputFolder, 'full_system_id_stribeck_result.mat');
parameterCsvFile = fullfile(outputFolder, 'full_system_id_stribeck_parameters.csv');
correlationCsvFile = fullfile(outputFolder, 'full_system_id_stribeck_correlations.csv');

%% ================================================================
% Start console logging
% ================================================================

if useDiary
    diary(fullfile(outputFolder, consoleLogFile));
end

fprintf('\n======================================================\n');
fprintf('Starting full-system Stribeck identification\n');
fprintf('Passive-link parameters are locked/fixed\n');
fprintf('======================================================\n\n');

fprintf('Identification runs: ');
fprintf('%d ', idxEst);
fprintf('\nValidation runs:     ');
fprintf('%d ', idxVal);
fprintf('\n\n');

modelPath = which(modelFile);
if isempty(modelPath)
    error(['MATLAB cannot find the model file "%s".\n' ...
           'Check modelsFolder and run: which %s'], modelFile, modelFile);
else
    fprintf('Using model file: %s\n\n', modelPath);
end

%% ================================================================
% Load locked passive-link parameters
% ================================================================

if loadPassiveParametersFromMat
    if ~isfile(passiveResultMatFile)
        error(['Passive-link result file not found:\n%s\n\n' ...
               'Update passiveResultMatFile near the top of this script, ' ...
               'or set loadPassiveParametersFromMat = false to use the manual fallback values.'], ...
               passiveResultMatFile);
    end

    fprintf('Loading passive-link parameters from:\n%s\n\n', passiveResultMatFile);

    passiveLoaded = load(passiveResultMatFile);
    passiveModel = [];

    if isfield(passiveLoaded, 'passiveID')
        if isfield(passiveLoaded.passiveID, 'mStrBest')
            passiveModel = passiveLoaded.passiveID.mStrBest;
        elseif isfield(passiveLoaded.passiveID, 'mCoul')
            passiveModel = passiveLoaded.passiveID.mCoul;
        elseif isfield(passiveLoaded.passiveID, 'mVisc')
            passiveModel = passiveLoaded.passiveID.mVisc;
        elseif isfield(passiveLoaded.passiveID, 'model_viscous')
            passiveModel = passiveLoaded.passiveID.model_viscous;
        end
    elseif isfield(passiveLoaded, 'mStrBest')
        passiveModel = passiveLoaded.mStrBest;
    elseif isfield(passiveLoaded, 'mCoul')
        passiveModel = passiveLoaded.mCoul;
    elseif isfield(passiveLoaded, 'mVisc')
        passiveModel = passiveLoaded.mVisc;
    elseif isfield(passiveLoaded, 'model_viscous')
        passiveModel = passiveLoaded.model_viscous;
    end

    if isempty(passiveModel)
        error(['Could not find a passive-link idnlgrey model in the result file.\n' ...
               'Expected passiveID.mStrBest, passiveID.mCoul, passiveID.mVisc, or passiveID.model_viscous.']);
    end

    p_b2_val = NaN;
    p_g2_val = NaN;
    p_c2_val = NaN;
    p_sdelta2_val = NaN;
    v_s2_val = NaN;
    eps_v2_val = NaN;

    for kPassivePar = 1:length(passiveModel.Parameters)
        passiveParName = char(passiveModel.Parameters(kPassivePar).Name);
        passiveParValue = passiveModel.Parameters(kPassivePar).Value;

        if iscell(passiveParValue)
            passiveParValue = passiveParValue{1};
        end

        passiveParValue = double(passiveParValue);

        switch passiveParName
            case 'p_b2'
                p_b2_val = passiveParValue;
            case 'p_g2'
                p_g2_val = passiveParValue;
            case 'p_c2'
                p_c2_val = passiveParValue;
            case {'p_sdelta2', 'p_sdelta'}
                p_sdelta2_val = passiveParValue;
            case {'v_s2', 'v_s'}
                v_s2_val = passiveParValue;
            case {'eps_v2', 'eps_v'}
                eps_v2_val = passiveParValue;
        end
    end

    if any(isnan([p_b2_val, p_g2_val, p_c2_val, p_sdelta2_val, v_s2_val, eps_v2_val]))
        error(['The selected passive-link result file does not contain all parameters ' ...
               'needed by the full Stribeck model.\n' ...
               'Required: p_b2, p_g2, p_c2, p_sdelta2, v_s2, eps_v2.\n' ...
               'Use the 3-step passive Stribeck result file, not the viscous-only passive result file.']);
    end
end

fprintf('Locked passive-link parameters used in the full-system model:\n');
fprintf('  p_b2      = %.9g\n', p_b2_val);
fprintf('  p_g2      = %.9g\n', p_g2_val);
fprintf('  p_c2      = %.9g\n', p_c2_val);
fprintf('  p_sdelta2 = %.9g\n', p_sdelta2_val);
fprintf('  v_s2      = %.9g\n', v_s2_val);
fprintf('  eps_v2    = %.9g\n\n', eps_v2_val);

%% ================================================================
% Load data and build multi-experiment iddata objects
% ================================================================

zEst = [];
zVal = [];

x0Est = [];
x0Val = [];

experimentNamesEst = {};
experimentNamesVal = {};

nEst = 0;
nVal = 0;

for kRun = 1:length(amplitudes)

    ampText = strrep(sprintf('%.2f', amplitudes(kRun)), '.', 'p');

    if runPrbs
        prbsFile = fullfile(prbsFolder, sprintf('fullsystem_prbs_A%s_run%02d.mat', ampText, kRun));

        if ~isfile(prbsFile)
            error('PRBS data file not found: %s', prbsFile);
        end

        load(prbsFile, 'theta_1', 'theta_2', 'u_ts');

        theta1 = double(squeeze(theta_1.Data(:)));
        theta2 = double(squeeze(theta_2.Data(:)));
        u = inputSign * double(squeeze(u_ts.Data(:)));

        y = [theta1, theta2];

        z = iddata(y, u, Ts);
        z.Name = sprintf('prbs_run_%02d', kRun);
        z.InputName = {'u'};
        z.InputUnit = {'V'};
        z.OutputName = {'theta_1', 'theta_2'};
        z.OutputUnit = {'rad', 'rad'};
        z.TimeUnit = 's';

        x0 = [theta1(1); theta2(1); 0; 0];

        if ismember(kRun, idxEst)
            nEst = nEst + 1;

            if isempty(zEst)
                zEst = z;
            else
                zEst = merge(zEst, z);
            end

            x0Est(:, nEst) = x0;
            experimentNamesEst{nEst, 1} = z.Name;
        else
            nVal = nVal + 1;

            if isempty(zVal)
                zVal = z;
            else
                zVal = merge(zVal, z);
            end

            x0Val(:, nVal) = x0;
            experimentNamesVal{nVal, 1} = z.Name;
        end
    end

    if runChirp
        chirpFile = fullfile(chirpFolder, sprintf('fullsystem_chirp_A%s_run%02d.mat', ampText, kRun));

        if ~isfile(chirpFile)
            error('Chirp data file not found: %s', chirpFile);
        end

        load(chirpFile, 'theta_1', 'theta_2', 'u_ts');

        theta1 = double(squeeze(theta_1.Data(:)));
        theta2 = double(squeeze(theta_2.Data(:)));
        u = inputSign * double(squeeze(u_ts.Data(:)));

        y = [theta1, theta2];

        z = iddata(y, u, Ts);
        z.Name = sprintf('chirp_run_%02d', kRun);
        z.InputName = {'u'};
        z.InputUnit = {'V'};
        z.OutputName = {'theta_1', 'theta_2'};
        z.OutputUnit = {'rad', 'rad'};
        z.TimeUnit = 's';

        x0 = [theta1(1); theta2(1); 0; 0];

        if ismember(kRun, idxEst)
            nEst = nEst + 1;

            if isempty(zEst)
                zEst = z;
            else
                zEst = merge(zEst, z);
            end

            x0Est(:, nEst) = x0;
            experimentNamesEst{nEst, 1} = z.Name;
        else
            nVal = nVal + 1;

            if isempty(zVal)
                zVal = z;
            else
                zVal = merge(zVal, z);
            end

            x0Val(:, nVal) = x0;
            experimentNamesVal{nVal, 1} = z.Name;
        end
    end
end

fprintf('Number of identification experiments: %d\n', nEst);
fprintf('Number of validation experiments:     %d\n\n', nVal);

%% ================================================================
% Grey-box model setup
% ================================================================

order = [2 1 4];     % 2 outputs, 1 input, 4 states

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

parameterNames = {'p_a', 'p_b1', 'p_c1', 'p_g1', 'p_u', 'p_0', ...
    'p_b2', 'p_g2', 'p_c2', 'p_sdelta2', 'v_s2', 'eps_v1', 'eps_v2'};

parameters = {
    p_a_init;
    p_b1_init;
    p_c1_init;
    p_g1_init;
    p_u_init;
    p_0_init;
    p_b2_val;
    p_g2_val;
    p_c2_val;
    p_sdelta2_val;
    v_s2_val;
    eps_v1_val;
    eps_v2_val
};

initialStatesEst = {
    x0Est(1, :);
    x0Est(2, :);
    x0Est(3, :);
    x0Est(4, :)
};

model0 = idnlgrey(modelFile, order, parameters, initialStatesEst, TsModel);

for kPar = 1:length(parameterNames)
    model0.Parameters(kPar).Name = parameterNames{kPar};
    model0.Parameters(kPar).Minimum = minimumValues(kPar);
    model0.Parameters(kPar).Maximum = maximumValues(kPar);
    model0.Parameters(kPar).Fixed = fixedStage1(kPar);
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

% Keep the initial states fixed during nlgreyest, as in the original script.
model0 = setinit(model0, 'Fixed', {
    true(1, nEst);
    true(1, nEst);
    true(1, nEst);
    true(1, nEst)
});

%% ================================================================
% Compare initial model before estimation
% ================================================================

if showInitialCompare
    figure('Name', 'Initial model fit before estimation', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, model0, compareOpt);
    title('Initial model fit before estimation');
end

%% ================================================================
% Stage 1: estimate active parameters with p_c1 fixed
% ================================================================

fprintf('\n======================================================\n');
fprintf('Stage 1: estimating active parameters with p_c1 fixed\n');
fprintf('Free parameters: p_a, p_b1, p_g1, p_u, p_0\n');
fprintf('Fixed parameters: p_c1 and all passive-link parameters\n');
fprintf('======================================================\n\n');

modelStage1 = nlgreyest(zEst, model0, optStage1);
modelStage1.Name = 'Stage 1 active model, p_c1 fixed';

fprintf('\nStage 1 parameter values:\n');
for kPar = 1:length(modelStage1.Parameters)
    if modelStage1.Parameters(kPar).Fixed
        fprintf('%-10s = %12.6g   fixed\n', modelStage1.Parameters(kPar).Name, modelStage1.Parameters(kPar).Value);
    else
        fprintf('%-10s = %12.6g   estimated\n', modelStage1.Parameters(kPar).Name, modelStage1.Parameters(kPar).Value);
    end
end

[~, fitStage1Est] = compare(zEst, modelStage1, compareOpt);

fitStage1Values = [];
if iscell(fitStage1Est)
    for kFit = 1:numel(fitStage1Est)
        fitStage1Values = [fitStage1Values; fitStage1Est{kFit}(:)]; %#ok<AGROW>
    end
else
    fitStage1Values = fitStage1Est(:);
end
fitStage1Values = fitStage1Values(isfinite(fitStage1Values));

fprintf('\nStage 1 mean identification fit: %.2f %%\n', mean(fitStage1Values));

if showStage1Compare
    figure('Name', 'Stage 1 identification fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, modelStage1, compareOpt);
    title('Stage 1 identification fit, p_{c1} fixed');
end

%% ================================================================
% Stage 2: unlock p_c1 and estimate final active parameter set
% ================================================================

fprintf('\n======================================================\n');
fprintf('Stage 2: estimating final active parameter set\n');
fprintf('Free parameters: p_a, p_b1, p_c1, p_g1, p_u, p_0\n');
fprintf('Fixed parameters: all passive-link parameters and smoothing constants\n');
fprintf('======================================================\n\n');

modelStage2Start = modelStage1;
modelStage2Start.Parameters(3).Value = stage2_p_c1_initial_value;

for kPar = 1:length(parameterNames)
    modelStage2Start.Parameters(kPar).Fixed = fixedStage2(kPar);
end

modelEst = nlgreyest(zEst, modelStage2Start, optStage2);
modelEst.Name = 'Final full-system Stribeck model with locked passive-link parameters';

%% ================================================================
% Validation model: same parameters, validation initial states
% ================================================================

initialStatesVal = {
    x0Val(1, :);
    x0Val(2, :);
    x0Val(3, :);
    x0Val(4, :)
};

modelVal = modelEst;
modelVal = setinit(modelVal, 'Value', initialStatesVal);
modelVal = setinit(modelVal, 'Fixed', {
    true(1, nVal);
    true(1, nVal);
    true(1, nVal);
    true(1, nVal)
});

%% ================================================================
% Compare final model on identification and validation data
% ================================================================

[yEst, fitEst] = compare(zEst, modelEst, compareOpt);
[yVal, fitVal] = compare(zVal, modelVal, compareOpt);

fitEstValues = [];
if iscell(fitEst)
    for kFit = 1:numel(fitEst)
        fitEstValues = [fitEstValues; fitEst{kFit}(:)]; %#ok<AGROW>
    end
else
    fitEstValues = fitEst(:);
end
fitEstValues = fitEstValues(isfinite(fitEstValues));

fitValValues = [];
if iscell(fitVal)
    for kFit = 1:numel(fitVal)
        fitValValues = [fitValValues; fitVal{kFit}(:)]; %#ok<AGROW>
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

if showFinalCompare
    figure('Name', 'Final identification fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, modelEst, compareOpt);
    title('Final identification fit');

    figure('Name', 'Final validation fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zVal, modelVal, compareOpt);
    title('Final validation fit');
end

%% ================================================================
% Parameter values, covariance, standard deviations and correlations
% ================================================================

fprintf('\n======================================================\n');
fprintf('Final parameter values and uncertainty\n');
fprintf('======================================================\n\n');

nPar = length(modelEst.Parameters);
parameter = cell(nPar, 1);
value = zeros(nPar, 1);
fixed = false(nPar, 1);
variance = NaN(nPar, 1);
stdDev = NaN(nPar, 1);
relativeStdPercent = NaN(nPar, 1);

for kPar = 1:nPar
    parameter{kPar} = modelEst.Parameters(kPar).Name;
    value(kPar) = modelEst.Parameters(kPar).Value;
    fixed(kPar) = modelEst.Parameters(kPar).Fixed;
end

freeIdx = find(~fixed);
covFree = getcov(modelEst, 'value', 'free');

for kFree = 1:length(freeIdx)
    variance(freeIdx(kFree)) = covFree(kFree, kFree);
    stdDev(freeIdx(kFree)) = sqrt(max(covFree(kFree, kFree), 0));

    if abs(value(freeIdx(kFree))) > eps
        relativeStdPercent(freeIdx(kFree)) = 100 * stdDev(freeIdx(kFree)) / abs(value(freeIdx(kFree)));
    end
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

fprintf('\nFinal loss function: %g\n', modelEst.Report.Fit.LossFcn);

%% ================================================================
% Save results
% ================================================================

save(resultMatFile, ...
    'modelEst', 'modelVal', 'modelStage1', 'zEst', 'zVal', ...
    'yEst', 'yVal', 'fitEst', 'fitVal', 'fitStage1Est', ...
    'resultTable', 'covFree', 'corrTable', ...
    'idxEst', 'idxVal', 'amplitudes', 'inputSign', ...
    'experimentNamesEst', 'experimentNamesVal', ...
    'loadPassiveParametersFromMat', 'passiveResultMatFile', ...
    'p_b2_val', 'p_g2_val', 'p_c2_val', 'p_sdelta2_val', 'v_s2_val', 'eps_v2_val');

writetable(resultTable, parameterCsvFile);
writetable(corrTable, correlationCsvFile, 'WriteRowNames', true);

fprintf('\nSaved result file: %s\n', resultMatFile);
fprintf('Saved parameter table: %s\n', parameterCsvFile);
fprintf('Saved correlation table: %s\n', correlationCsvFile);

if useDiary
    diary off;
end
