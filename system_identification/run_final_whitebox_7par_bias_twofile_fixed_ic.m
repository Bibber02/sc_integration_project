clear;
clear functions;
clc;
close all;

fprintf('Running run_final_whitebox_7par_bias_twofile_fixed_ic.m\n');

%% ================================================================
% User settings
% ================================================================

% These may be .mat files containing theta_1, theta_2, u_ts,
% OR .m scripts that create theta_1, theta_2, u_ts in the workspace.trainFile = "identification_data\whole_system_identification\prbs_train_121s_amp_015.mat";
valFile   = "identification_data\whole_system_identification\prbs_validation_121s_amp_015.mat";

% Variables inside each file/script.
% theta_1 and theta_2 are expected to be timeseries.
% u_ts is expected to be timeseries.
theta1VarName = "theta_1";
theta2VarName = "theta_2";
uVarName      = "u_ts";

% The data has 121 seconds total with the first second idle.
usedEndTime = 121.0;
idleStart   = 0.0;
idleEnd     = 1.0;

% Sign conventions from earlier tests.
inputPolarity = -1;
relativeAngleSign = 1;

% Since the first second is idle, subtract the mean over 0-1 s.
% This is usually better than subtracting only the first sample.
removeIdleMeanOffsets = true;
removeIdleMeanInput   = false;   % leave false unless u_ts has nonzero idle bias

% Fixed initial conditions.
% Because the first second is idle and offsets are removed using the idle mean,
% use x0 = [0; 0; 0; 0].
useFixedZeroInitialConditions = false;

% Output weighting during parameter estimation.
outputWeight = diag([1, 8]);

% Estimation settings.
parameterMaxIterations = 150;

%% ================================================================
% Load training and validation runs
% ================================================================

trainRun = loadOneRun(trainFile, theta1VarName, theta2VarName, uVarName);
valRun   = loadOneRun(valFile,   theta1VarName, theta2VarName, uVarName);

trainRun = preprocessRun(trainRun, usedEndTime, idleStart, idleEnd, ...
    inputPolarity, relativeAngleSign, removeIdleMeanOffsets, removeIdleMeanInput, "training");

valRun = preprocessRun(valRun, usedEndTime, idleStart, idleEnd, ...
    inputPolarity, relativeAngleSign, removeIdleMeanOffsets, removeIdleMeanInput, "validation");

% Check sample times.
Ts = trainRun.Ts;

if abs(valRun.Ts - Ts) > 1e-9
    warning('Training Ts and validation Ts differ: train %.12g, val %.12g', Ts, valRun.Ts);
end

fprintf('\nTraining samples:   %d, duration %.3f s\n', numel(trainRun.u), trainRun.t(end));
fprintf('Validation samples: %d, duration %.3f s\n', numel(valRun.u), valRun.t(end));
fprintf('Sampling time Ts: %.8f s\n', Ts);
fprintf('Sampling frequency fs: %.3f Hz\n', 1/Ts);

plotMeasuredRuns(trainRun, valRun);

%% ================================================================
% Build iddata objects
% ================================================================

zTrain = iddata([trainRun.theta1, trainRun.phi2], trainRun.u, Ts);
zVal   = iddata([valRun.theta1,   valRun.phi2],   valRun.u,   Ts);

zTrain.InputName  = {'u'};
zTrain.OutputName = {'theta_1', 'phi_2'};
zTrain.InputUnit  = {'input'};
zTrain.OutputUnit = {'rad', 'rad'};
zTrain.TimeUnit   = 's';

zVal.InputName  = {'u'};
zVal.OutputName = {'theta_1', 'phi_2'};
zVal.InputUnit  = {'input'};
zVal.OutputUnit = {'rad', 'rad'};
zVal.TimeUnit   = 's';

%% ================================================================
% Create model
% ================================================================

modelFile = 'final_whitebox_7par_bias_fixed_ic_twofile_model';

Order = [2 1 4];
TsModel = 0;

% Initial parameter guesses from previous successful 6/7 parameter fits.
p0.p_J1    = 14.6;
p0.p_kappa = 0.32;
p0.p_g1    = 97.0;
p0.p_g2    = 112.0;
p0.p_Ku    = 4243.0;
p0.p_b1    = 756.0;
p0.p_tau0  = -4.25;

if useFixedZeroInitialConditions
    InitialStates = {0; 0; 0; 0};
else
    theta1_0 = trainRun.theta1(1);
    phi2_0   = trainRun.phi2(1);
    theta2_0 = theta1_0 + phi2_0;
    InitialStates = {theta1_0; theta2_0; 0; 0};
end

Parameters = {
    p0.p_J1
    p0.p_kappa
    p0.p_g1
    p0.p_g2
    p0.p_Ku
    p0.p_b1
    p0.p_tau0
};

model0 = idnlgrey(modelFile, Order, Parameters, InitialStates, TsModel);
model0 = configureModel(model0);

fprintf('\nInitial model:\n');
dispParameters(model0);
dispInitialStates(model0);

%% ================================================================
% Estimate parameters on training run
% ================================================================

opt = nlgreyestOptions;
opt.Display = 'on';
opt.EstimateCovariance = true;
opt.SearchOptions.MaxIterations = parameterMaxIterations;
opt.OutputWeight = outputWeight;

fprintf('\n============================================================\n');
fprintf('Estimate 7 parameters on training run, fixed initial conditions\n');
fprintf('============================================================\n');

modelEst = nlgreyest(zTrain, model0, opt);
modelEst.Name = 'Final 7-parameter bias model, two-file fixed-IC fit';

fprintf('\nEstimated parameters:\n');
dispParameters(modelEst);

fprintf('Initial states used during parameter estimation:\n');
dispInitialStates(modelEst);

%% ================================================================
% Compare train and validation with fixed zero initial conditions
% ================================================================

modelTrain = setFixedInitialStateForRun(modelEst, trainRun);
modelVal   = setFixedInitialStateForRun(modelEst, valRun);

[yHatTrain, fitTrain] = compareUsingModelIC(zTrain, modelTrain);
[yHatVal,   fitVal]   = compareUsingModelIC(zVal,   modelVal);

fprintf('\n============================================================\n');
fprintf('FINAL FIT WITH FIXED INITIAL CONDITIONS\n');
fprintf('============================================================\n');

fprintf('\nTraining fit:\n');
dispFitBreakdown(fitTrain, {'theta_1', 'phi_2'});

fprintf('Validation fit:\n');
dispFitBreakdown(fitVal, {'theta_1', 'phi_2'});

fprintf('\nFinal parameter values:\n');
dispParameters(modelEst);

plotFinalComparison(zTrain, zVal, yHatTrain, yHatVal, trainRun.t, valRun.t, fitTrain, fitVal);

fprintf('\nAmplitude diagnostics, training:\n');
printAmplitudeDiagnostics(zTrain, yHatTrain, {'theta_1', 'phi_2'});

fprintf('\nAmplitude diagnostics, validation:\n');
printAmplitudeDiagnostics(zVal, yHatVal, {'theta_1', 'phi_2'});

%% ================================================================
% Save result
% ================================================================

finalWhitebox7FixedIC.modelEst = modelEst;
finalWhitebox7FixedIC.modelTrain = modelTrain;
finalWhitebox7FixedIC.modelVal = modelVal;
finalWhitebox7FixedIC.fitTrain = fitTrain;
finalWhitebox7FixedIC.fitVal = fitVal;
finalWhitebox7FixedIC.trainRun = trainRun;
finalWhitebox7FixedIC.valRun = valRun;

finalWhitebox7FixedIC.settings.trainFile = trainFile;
finalWhitebox7FixedIC.settings.valFile = valFile;
finalWhitebox7FixedIC.settings.theta1VarName = theta1VarName;
finalWhitebox7FixedIC.settings.theta2VarName = theta2VarName;
finalWhitebox7FixedIC.settings.uVarName = uVarName;
finalWhitebox7FixedIC.settings.inputPolarity = inputPolarity;
finalWhitebox7FixedIC.settings.relativeAngleSign = relativeAngleSign;
finalWhitebox7FixedIC.settings.removeIdleMeanOffsets = removeIdleMeanOffsets;
finalWhitebox7FixedIC.settings.removeIdleMeanInput = removeIdleMeanInput;
finalWhitebox7FixedIC.settings.outputWeight = outputWeight;

save('final_whitebox_7par_bias_twofile_fixed_ic_result.mat', 'finalWhitebox7FixedIC');

fprintf('\nSaved result to final_whitebox_7par_bias_twofile_fixed_ic_result.mat\n');

%% ================================================================
% Local helper functions
% ================================================================

function runData = loadOneRun(fileName, theta1VarName, theta2VarName, uVarName)

fprintf('\nLoading run file: %s\n', fileName);

fileNameChar = char(fileName);
[~, ~, ext] = fileparts(fileNameChar);

switch lower(ext)
    case '.mat'
        S = load(fileNameChar);

    case '.m'
        % Run script in an isolated function workspace and collect the
        % required variables afterwards.
        S = struct();
        run(fileNameChar);

        if exist(theta1VarName, 'var')
            S.(theta1VarName) = eval(theta1VarName);
        end

        if exist(theta2VarName, 'var')
            S.(theta2VarName) = eval(theta2VarName);
        end

        if exist(uVarName, 'var')
            S.(uVarName) = eval(uVarName);
        end

    otherwise
        error('Unsupported file type: %s. Use .mat or .m.', ext);
end

required = [theta1VarName, theta2VarName, uVarName];

for k = 1:numel(required)
    if ~isfield(S, required(k))
        error('File %s does not contain required variable "%s".', fileName, required(k));
    end
end

theta1_ts = S.(theta1VarName);
theta2_ts = S.(theta2VarName);
u_ts      = S.(uVarName);

runData.t_raw      = double(theta1_ts.Time(:));
runData.theta1_raw = double(theta1_ts.Data(:));
runData.theta2_raw = double(theta2_ts.Data(:));
runData.u_raw      = double(u_ts.Data(:));

N = min([numel(runData.t_raw), numel(runData.theta1_raw), ...
         numel(runData.theta2_raw), numel(runData.u_raw)]);

runData.t_raw      = runData.t_raw(1:N);
runData.theta1_raw = runData.theta1_raw(1:N);
runData.theta2_raw = runData.theta2_raw(1:N);
runData.u_raw      = runData.u_raw(1:N);

runData.t_raw = runData.t_raw - runData.t_raw(1);

end

function runData = preprocessRun(runData, usedEndTime, idleStart, idleEnd, ...
    inputPolarity, relativeAngleSign, removeIdleMeanOffsets, removeIdleMeanInput, label)

idxUsed = runData.t_raw <= usedEndTime;

t      = runData.t_raw(idxUsed);
theta1 = runData.theta1_raw(idxUsed);
phi2   = relativeAngleSign * runData.theta2_raw(idxUsed);
u      = inputPolarity * runData.u_raw(idxUsed);

t = t - t(1);

idxIdle = t >= idleStart & t <= idleEnd;

if ~any(idxIdle)
    warning('No idle samples found for %s run. Using first sample for offsets.', label);
    idxIdle = false(size(t));
    idxIdle(1) = true;
end

if removeIdleMeanOffsets
    theta1Offset = mean(theta1(idxIdle), 'omitnan');
    phi2Offset   = mean(phi2(idxIdle),   'omitnan');

    theta1 = theta1 - theta1Offset;
    phi2   = phi2   - phi2Offset;
else
    theta1Offset = 0;
    phi2Offset   = 0;
end

if removeIdleMeanInput
    uOffset = mean(u(idxIdle), 'omitnan');
    u = u - uOffset;
else
    uOffset = 0;
end

runData.t = t;
runData.theta1 = theta1;
runData.phi2 = phi2;
runData.theta2_abs = theta1 + phi2;
runData.u = u;

runData.Ts = mean(diff(t));
runData.fs = 1 / runData.Ts;

runData.theta1Offset = theta1Offset;
runData.phi2Offset = phi2Offset;
runData.uOffset = uOffset;
runData.label = label;

fprintf('\n%s run preprocessing:\n', label);
fprintf('  duration: %.3f s\n', t(end));
fprintf('  samples:  %d\n', numel(t));
fprintf('  idle window: %.3f to %.3f s\n', idleStart, idleEnd);
fprintf('  theta_1 idle offset removed: %.8g rad\n', theta1Offset);
fprintf('  phi_2 idle offset removed:   %.8g rad\n', phi2Offset);
fprintf('  input idle offset removed:   %.8g\n', uOffset);
fprintf('  first processed theta_1: %.8g rad\n', theta1(1));
fprintf('  first processed phi_2:   %.8g rad\n', phi2(1));

end

function model = configureModel(model)

names = {
    'p_J1'
    'p_kappa'
    'p_g1'
    'p_g2'
    'p_Ku'
    'p_b1'
    'p_tau0'
};

units = {
    '-'
    '-'
    '1/s^2'
    '1/s^2'
    'torque/input'
    '1/s'
    'torque'
};

mins = [
    0.05
   -0.98
    0
    0
    0
    0
   -5000
];

maxs = [
    80
    0.98
    800
    800
    20000
    5000
    5000
];

for k = 1:numel(names)
    model.Parameters(k).Name = names{k};
    model.Parameters(k).Unit = units{k};
    model.Parameters(k).Minimum = mins(k);
    model.Parameters(k).Maximum = maxs(k);
    model.Parameters(k).Fixed = false;
end

model.TimeUnit = 's';
model.InputName = {'u'};
model.InputUnit = {'input'};
model.OutputName = {'theta_1', 'phi_2'};
model.OutputUnit = {'rad', 'rad'};

model.InitialStates(1).Name = 'theta_1';
model.InitialStates(1).Unit = 'rad';
model.InitialStates(2).Name = 'theta_2_abs';
model.InitialStates(2).Unit = 'rad';
model.InitialStates(3).Name = 'omega_1';
model.InitialStates(3).Unit = 'rad/s';
model.InitialStates(4).Name = 'omega_2_abs';
model.InitialStates(4).Unit = 'rad/s';

% Fixed initial conditions: no IC estimation.
model = setinit(model, 'Value', {0; 0; 0; 0});
model = setinit(model, 'Fixed', {true; true; true; true});

end

function modelOut = setFixedInitialStateForRun(modelIn, runData)

modelOut = modelIn;

% Because each run has an idle first second and mean offsets are removed,
% the physical fixed initial state is zero.
modelOut = setinit(modelOut, 'Value', {0; 0; 0; 0});
modelOut = setinit(modelOut, 'Fixed', {true; true; true; true});

end

function [yHat, fit] = compareUsingModelIC(data, model)

try
    optCompare = compareOptions;
    optCompare.InitialCondition = 'model';
    [yHat, fit] = compare(data, model, optCompare);
catch
    try
        [yHat, fit] = compare(data, model, 'InitialCondition', 'model');
    catch
        [yHat, fit] = compare(data, model);
    end
end

end

function plotMeasuredRuns(trainRun, valRun)

figure('Name', 'Measured training and validation runs', ...
    'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.8]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(trainRun.t, trainRun.theta1);
grid on;
ylabel('\theta_1 [rad]');
title('Training \theta_1');

nexttile;
plot(valRun.t, valRun.theta1);
grid on;
ylabel('\theta_1 [rad]');
title('Validation \theta_1');

nexttile;
plot(trainRun.t, trainRun.phi2);
grid on;
ylabel('\phi_2 [rad]');
title('Training \phi_2');

nexttile;
plot(valRun.t, valRun.phi2);
grid on;
ylabel('\phi_2 [rad]');
title('Validation \phi_2');

nexttile;
plot(trainRun.t, trainRun.u);
grid on;
ylabel('u');
xlabel('Time [s]');
title('Training input');

nexttile;
plot(valRun.t, valRun.u);
grid on;
ylabel('u');
xlabel('Time [s]');
title('Validation input');

end

function plotFinalComparison(zTrain, zVal, yHatTrain, yHatVal, tTrain, tVal, fitTrain, fitVal)

yTrain = getOutputDataMatrix(zTrain);
yVal   = getOutputDataMatrix(zVal);

yTrainHat = getOutputDataMatrix(yHatTrain);
yValHat   = getOutputDataMatrix(yHatVal);

trainFitMean = meanFitValue(fitTrain);
valFitMean   = meanFitValue(fitVal);

figure('Name', 'Final 7-parameter fixed-IC train/validation fit', ...
    'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(tTrain, yTrain(:,1), 'DisplayName', 'measured');
hold on;
nTrain1 = min(numel(tTrain), size(yTrainHat, 1));
plot(tTrain(1:nTrain1), yTrainHat(1:nTrain1,1), '--', 'DisplayName', 'model');
grid on;
ylabel('\theta_1 [rad]');
title(sprintf('Training \\theta_1, mean %.2f %%', trainFitMean));
legend('Location', 'best');

nexttile;
plot(tTrain, yTrain(:,2), 'DisplayName', 'measured');
hold on;
nTrain2 = min(numel(tTrain), size(yTrainHat, 1));
plot(tTrain(1:nTrain2), yTrainHat(1:nTrain2,2), '--', 'DisplayName', 'model');
grid on;
ylabel('\phi_2 [rad]');
title('Training \phi_2');
legend('Location', 'best');

nexttile;
plot(tVal, yVal(:,1), 'DisplayName', 'measured');
hold on;
nVal1 = min(numel(tVal), size(yValHat, 1));
plot(tVal(1:nVal1), yValHat(1:nVal1,1), '--', 'DisplayName', 'model');
grid on;
ylabel('\theta_1 [rad]');
xlabel('Time [s]');
title(sprintf('Validation \\theta_1, mean %.2f %%', valFitMean));
legend('Location', 'best');

nexttile;
plot(tVal, yVal(:,2), 'DisplayName', 'measured');
hold on;
nVal2 = min(numel(tVal), size(yValHat, 1));
plot(tVal(1:nVal2), yValHat(1:nVal2,2), '--', 'DisplayName', 'model');
grid on;
ylabel('\phi_2 [rad]');
xlabel('Time [s]');
title('Validation \phi_2');
legend('Location', 'best');

sgtitle('Final 7-parameter bias model, two separate runs, fixed ICs');

end

function y = getOutputDataMatrix(dataObject)

try
    y = get(dataObject, 'OutputData');
catch
    y = dataObject.OutputData;
end

if iscell(y)
    y = y{1};
end

y = double(y);

if isvector(y)
    y = y(:);
end

end

function dispParameters(model)

for k = 1:length(model.Parameters)
    fprintf('  %-10s = %.8g', model.Parameters(k).Name, model.Parameters(k).Value);

    if isParameterFixed(model.Parameters(k))
        fprintf('  fixed\n');
    else
        fprintf('  free\n');
    end
end

end

function dispInitialStates(model)

for k = 1:length(model.InitialStates)
    fprintf('  %-12s = %.8g', model.InitialStates(k).Name, model.InitialStates(k).Value);

    fixed = model.InitialStates(k).Fixed;
    if iscell(fixed)
        fixed = fixed{1};
    end

    if all(fixed(:))
        fprintf('  fixed\n');
    else
        fprintf('  free\n');
    end
end

end

function tf = isParameterFixed(parameter)

tf = parameter.Fixed;

if iscell(tf)
    tf = tf{1};
end

tf = all(tf(:));

end

function dispFitBreakdown(fit, outputNames)

values = collectFitValues(fit);
values = values(isfinite(values));

if isempty(values)
    fprintf('  mean      = NaN %%\n');
    return;
end

for k = 1:min(numel(values), numel(outputNames))
    fprintf('  %-8s = %.3f %%\n', outputNames{k}, values(k));
end

fprintf('  %-8s = %.3f %%\n', 'mean', mean(values));

end

function m = meanFitValue(fit)

values = collectFitValues(fit);
values = values(isfinite(values));

if isempty(values)
    m = NaN;
else
    m = mean(values);
end

end

function values = collectFitValues(fit)

if isnumeric(fit)
    values = fit(:);
elseif iscell(fit)
    values = [];

    for k = 1:numel(fit)
        values = [values; collectFitValues(fit{k})]; %#ok<AGROW>
    end
else
    values = NaN;
end

end

function printAmplitudeDiagnostics(dataMeasured, dataModel, outputNames)

ym = getOutputDataMatrix(dataMeasured);
yh = getOutputDataMatrix(dataModel);

nOut = min(size(ym,2), size(yh,2));

for k = 1:nOut
    measuredP2P = max(ym(:,k)) - min(ym(:,k));
    modelP2P    = max(yh(:,k)) - min(yh(:,k));

    if measuredP2P > 1e-12
        ratio = modelP2P / measuredP2P;
    else
        ratio = NaN;
    end

    fprintf('  %-8s measured p2p = %.6g, model p2p = %.6g, ratio = %.4f\n', ...
        outputNames{k}, measuredP2P, modelP2P, ratio);
end

end
