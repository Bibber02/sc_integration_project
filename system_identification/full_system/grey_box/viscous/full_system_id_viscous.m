clear;
clear functions;
clc;
close all;

%% ================================================================
% Full-system nonlinear grey-box identification, viscous-only version
%
% This version keeps the same data handling and validation structure as the
% Stribeck full-system script, but removes all Coulomb/Stribeck friction
% terms from the full-system model.
%
% The model only contains:
%   Joint 1 friction: p_b1 * theta_1_dot
%   Joint 2 friction: p_b2 * theta_2_dot
%
% The script first estimates the active/full-system parameters while keeping
% the passive-link parameters p_b2 and p_g2 fixed to the passive-link ID result.
% It then performs a second estimation round where p_b2 and p_g2 are unlocked.
% This gives the full-system data one additional chance to refine the passive
% viscous parameters without changing the initial locked-passive baseline.
%
% Required model file on the MATLAB path:
%   greybox_id_full_viscous_model.m
% ================================================================

%% ================================================================
% Configuration
% ================================================================

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = scriptFolder;
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(projectRoot);
scip.setupPath;
projectPaths = scip.paths;

% Change these paths if your folder structure is different.
dataFolder   = projectPaths.fullSystemMeasurementData;
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

modelFile = 'greybox_id_full_viscous_model';

Ts = 0.01;
TsModel = 0;         % continuous-time grey-box model
inputSign = -1;      % original script used u = -u. Try +1 if the fit is mirrored.

amplitudes = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];

idxEst = [1 3 5 7 9];    % odd runs for identification
idxVal = [2 4 6 8 10];   % even runs for validation

runPrbs = true;
runChirp = true;

useDiary = true;
consoleLogFile = 'full_system_id_viscous_console_output.txt';

showInitialCompare = true;
showFinalCompare = true;

figureUnits = 'pixels';
figurePosition = [80 80 1150 720];

% Passive-link parameters are loaded from the passive-link ID result file.
% The viscous full-system model only needs p_b2 and p_g2.
% Keep the manual values as a fallback only.
loadPassiveParametersFromMat = true;
passiveResultMatFile = projectPaths.passiveLinkViscousResult;

% Manual fallback values. These are only used when
% loadPassiveParametersFromMat = false.
p_b2_val_manual = 0.067028846;
p_g2_val_manual = 112.02675;

% These variables are assigned after the model-file check below.
p_b2_val = p_b2_val_manual;
p_g2_val = p_g2_val_manual;

% Initial full-system parameter guesses.
p_a_init  = 14;
p_b1_init = 670;
p_g1_init = 97;
p_u_init  = 4000;
p_0_init  = 0;

% Bounds for the viscous full-system model parameters:
% [p_a, p_b1, p_g1, p_u, p_0, p_b2, p_g2]
minimumValues = [5, 100, 5, 500, -200, 0, 0];
maximumValues = [80, 4000, 500, 30000, 200, Inf, Inf];

% Round 1: estimate active/full-system parameters and keep passive-link parameters fixed.
fixedParametersRound1 = [false, false, false, false, false, true, true];

% Round 2: start from Round 1, then unlock p_b2 and p_g2 as an additional refinement step.
runSecondRoundLoosenPassive = true;
fixedParametersRound2 = [false, false, false, false, false, false, false];

% Estimation options.
optRound1 = nlgreyestOptions;
optRound1.Display = 'Full';
optRound1.EstimateCovariance = true;
optRound1.SearchMethod = 'lm';
optRound1.SearchOptions.MaxIterations = 150;
optRound1.OutputWeight = diag([1, 8]);

optRound2 = nlgreyestOptions;
optRound2.Display = 'Full';
optRound2.EstimateCovariance = true;
optRound2.SearchMethod = 'lm';
optRound2.SearchOptions.MaxIterations = 100;
optRound2.OutputWeight = diag([1, 8]);

compareOpt = compareOptions;
compareOpt.InitialCondition = 'estimate';

resultMatFile = fullfile(outputFolder, 'full_system_id_viscous_result.mat');
parameterCsvFile = fullfile(outputFolder, 'full_system_id_viscous_parameters.csv');
correlationCsvFile = fullfile(outputFolder, 'full_system_id_viscous_correlations.csv');
lockedParameterCsvFile = fullfile(outputFolder, 'full_system_id_viscous_parameters.csv');

%% ================================================================
% Start console logging
% ================================================================

if useDiary
    diary(fullfile(outputFolder, consoleLogFile));
end

fprintf('\n======================================================\n');
fprintf('Starting full-system viscous-only identification\n');
fprintf('Round 1 locks passive-link p_b2 and p_g2\n');
fprintf('Round 2 optionally unlocks passive-link p_b2 and p_g2\n');
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
        if isfield(passiveLoaded.passiveID, 'model_viscous')
            passiveModel = passiveLoaded.passiveID.model_viscous;
        elseif isfield(passiveLoaded.passiveID, 'mVisc')
            passiveModel = passiveLoaded.passiveID.mVisc;
        elseif isfield(passiveLoaded.passiveID, 'mStrBest')
            passiveModel = passiveLoaded.passiveID.mStrBest;
        elseif isfield(passiveLoaded.passiveID, 'mCoul')
            passiveModel = passiveLoaded.passiveID.mCoul;
        end
    elseif isfield(passiveLoaded, 'model_viscous')
        passiveModel = passiveLoaded.model_viscous;
    elseif isfield(passiveLoaded, 'mVisc')
        passiveModel = passiveLoaded.mVisc;
    elseif isfield(passiveLoaded, 'mStrBest')
        passiveModel = passiveLoaded.mStrBest;
    elseif isfield(passiveLoaded, 'mCoul')
        passiveModel = passiveLoaded.mCoul;
    end

    if isempty(passiveModel)
        error(['Could not find a passive-link idnlgrey model in the result file.\n' ...
               'Expected passiveID.model_viscous, passiveID.mVisc, passiveID.mStrBest, or passiveID.mCoul.']);
    end

    p_b2_val = NaN;
    p_g2_val = NaN;

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
        end
    end

    if any(isnan([p_b2_val, p_g2_val]))
        error(['The selected passive-link result file does not contain the parameters ' ...
               'needed by the full viscous model. Required: p_b2 and p_g2.']);
    end
end

fprintf('Locked passive-link parameters used in the full-system viscous model:\n');
fprintf('  p_b2 = %.9g\n', p_b2_val);
fprintf('  p_g2 = %.9g\n\n', p_g2_val);

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
%  1 p_a    active/full-system inertia ratio
%  2 p_b1   active joint effective viscous damping
%  3 p_g1   active-link gravity parameter
%  4 p_u    motor input gain
%  5 p_0    constant torque/input bias
%  6 p_b2   passive joint viscous damping, fixed
%  7 p_g2   passive-link gravity parameter, fixed

parameterNames = {'p_a', 'p_b1', 'p_g1', 'p_u', 'p_0', 'p_b2', 'p_g2'};

parameters = {
    p_a_init;
    p_b1_init;
    p_g1_init;
    p_u_init;
    p_0_init;
    p_b2_val;
    p_g2_val
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
    model0.Parameters(kPar).Fixed = fixedParametersRound1(kPar);
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
    figure('Name', 'Initial viscous model fit before estimation', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, model0, compareOpt);
    title('Initial viscous model fit before estimation');
end

%% ================================================================
% Round 1: estimate with passive-link parameters locked
% ================================================================

fprintf('\n======================================================\n');
fprintf('Round 1: viscous-only full-system model, locked passive-link parameters\n');
fprintf('Free parameters: p_a, p_b1, p_g1, p_u, p_0\n');
fprintf('Fixed parameters: p_b2, p_g2\n');
fprintf('======================================================\n\n');

modelLocked = nlgreyest(zEst, model0, optRound1);
modelLocked.Name = 'Round 1 full-system viscous-only model, locked passive-link parameters';

fprintf('\nRound 1 parameter values:\n');
for kPar = 1:length(modelLocked.Parameters)
    if modelLocked.Parameters(kPar).Fixed
        fprintf('%-10s = %12.6g   fixed\n', modelLocked.Parameters(kPar).Name, modelLocked.Parameters(kPar).Value);
    else
        fprintf('%-10s = %12.6g   estimated\n', modelLocked.Parameters(kPar).Name, modelLocked.Parameters(kPar).Value);
    end
end

[yLockedEst, fitLockedEst] = compare(zEst, modelLocked, compareOpt);

fitLockedEstValues = [];
if iscell(fitLockedEst)
    for kFit = 1:numel(fitLockedEst)
        fitLockedEstValues = [fitLockedEstValues; fitLockedEst{kFit}(:)]; %#ok<AGROW>
    end
else
    fitLockedEstValues = fitLockedEst(:);
end
fitLockedEstValues = fitLockedEstValues(isfinite(fitLockedEstValues));

fprintf('\nRound 1 mean identification fit: %.2f %%\n', mean(fitLockedEstValues));

if showFinalCompare
    figure('Name', 'Round 1 viscous identification fit, locked passive', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, modelLocked, compareOpt);
    title('Round 1 viscous identification fit, locked passive parameters');
end

%% ================================================================
% Round 2: unlock passive-link viscous parameters and estimate again
% ================================================================

if runSecondRoundLoosenPassive
    fprintf('\n======================================================\n');
    fprintf('Round 2: viscous-only full-system model, unlocked passive-link parameters\n');
    fprintf('Free parameters: p_a, p_b1, p_g1, p_u, p_0, p_b2, p_g2\n');
    fprintf('Starting point: Round 1 locked-passive model\n');
    fprintf('======================================================\n\n');

    modelUnlockedStart = modelLocked;

    for kPar = 1:length(parameterNames)
        modelUnlockedStart.Parameters(kPar).Fixed = fixedParametersRound2(kPar);
    end

    modelEst = nlgreyest(zEst, modelUnlockedStart, optRound2);
    modelEst.Name = 'Round 2 full-system viscous-only model, passive parameters unlocked';
else
    fprintf('\nSkipping Round 2 because runSecondRoundLoosenPassive = false.\n');
    modelUnlockedStart = modelLocked;
    modelEst = modelLocked;
    modelEst.Name = 'Final full-system viscous-only model, locked passive-link parameters';
end

fprintf('\nFinal parameter values after selected estimation round:\n');
for kPar = 1:length(modelEst.Parameters)
    if modelEst.Parameters(kPar).Fixed
        fprintf('%-10s = %12.6g   fixed\n', modelEst.Parameters(kPar).Name, modelEst.Parameters(kPar).Value);
    else
        fprintf('%-10s = %12.6g   estimated\n', modelEst.Parameters(kPar).Name, modelEst.Parameters(kPar).Value);
    end
end

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

modelLockedVal = modelLocked;
modelLockedVal = setinit(modelLockedVal, 'Value', initialStatesVal);
modelLockedVal = setinit(modelLockedVal, 'Fixed', {
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
[yLockedVal, fitLockedVal] = compare(zVal, modelLockedVal, compareOpt);

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

fitLockedValValues = [];
if iscell(fitLockedVal)
    for kFit = 1:numel(fitLockedVal)
        fitLockedValValues = [fitLockedValValues; fitLockedVal{kFit}(:)]; %#ok<AGROW>
    end
else
    fitLockedValValues = fitLockedVal(:);
end
fitLockedValValues = fitLockedValValues(isfinite(fitLockedValValues));

fprintf('\nRound 1 locked-passive validation fit, per experiment and output:\n');
disp(fitLockedVal);
fprintf('\nRound 1 mean validation fit: %.2f %%\n', mean(fitLockedValValues));

fprintf('\nFinal identification fit, per experiment and output:\n');
disp(fitEst);

fprintf('\nFinal validation fit, per experiment and output:\n');
disp(fitVal);

fprintf('\nFinal mean identification fit: %.2f %%\n', mean(fitEstValues));
fprintf('Final mean validation fit:     %.2f %%\n', mean(fitValValues));

if showFinalCompare
    figure('Name', 'Round 1 viscous validation fit, locked passive', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zVal, modelLockedVal, compareOpt);
    title('Round 1 viscous validation fit, locked passive parameters');

    figure('Name', 'Final viscous identification fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, modelEst, compareOpt);
    title('Final viscous identification fit');

    figure('Name', 'Final viscous validation fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zVal, modelVal, compareOpt);
    title('Final viscous validation fit');
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

nParLocked = length(modelLocked.Parameters);
lockedParameter = cell(nParLocked, 1);
lockedValue = zeros(nParLocked, 1);
lockedFixed = false(nParLocked, 1);

for kPar = 1:nParLocked
    lockedParameter{kPar} = modelLocked.Parameters(kPar).Name;
    lockedValue(kPar) = modelLocked.Parameters(kPar).Value;
    lockedFixed(kPar) = modelLocked.Parameters(kPar).Fixed;
end

lockedResultTable = table(lockedParameter, lockedValue, lockedFixed, ...
    'VariableNames', {'Parameter', 'Value', 'Fixed'});

fprintf('\nRound 1 locked-passive parameter table:\n');
disp(lockedResultTable);

fprintf('\nRound 1 locked-passive loss function: %g\n', modelLocked.Report.Fit.LossFcn);
fprintf('Final loss function: %g\n', modelEst.Report.Fit.LossFcn);

%% ================================================================
% Save results
% ================================================================

save(resultMatFile, ...
    'modelEst', 'modelVal', 'modelLocked', 'modelLockedVal', 'modelUnlockedStart', ...
    'zEst', 'zVal', ...
    'yEst', 'yVal', 'yLockedEst', 'yLockedVal', ...
    'fitEst', 'fitVal', 'fitLockedEst', 'fitLockedVal', ...
    'resultTable', 'lockedResultTable', 'covFree', 'corrTable', ...
    'idxEst', 'idxVal', 'amplitudes', 'inputSign', ...
    'experimentNamesEst', 'experimentNamesVal', ...
    'runSecondRoundLoosenPassive', ...
    'loadPassiveParametersFromMat', 'passiveResultMatFile', ...
    'p_b2_val', 'p_g2_val');

writetable(resultTable, parameterCsvFile);
writetable(lockedResultTable, lockedParameterCsvFile);
writetable(corrTable, correlationCsvFile, 'WriteRowNames', true);

fprintf('\nSaved result file: %s\n', resultMatFile);
fprintf('Saved final parameter table: %s\n', parameterCsvFile);
fprintf('Saved Round 1 locked-passive parameter table: %s\n', lockedParameterCsvFile);
fprintf('Saved correlation table: %s\n', correlationCsvFile);

if useDiary
    diary off;
end
