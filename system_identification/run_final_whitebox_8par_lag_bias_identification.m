clear;
clear functions;
clc;
close all;

fprintf('Running run_final_whitebox_8par_lag_bias_identification.m\n');

%% ================================================================
% User settings
% ================================================================

dataFile = "identification_data\whole_system_identification\prbs_long.mat";
datasetName = "PRBS_Long2_Amp0_15";

usedEndTime = 90.0;

trainStart = 0.0;
trainEnd   = 60.0;
valStart   = 60.0;
valEnd     = 90.0;

% Sign conventions from earlier tests.
inputPolarity = -1;
relativeAngleSign = 1;

% Treat the first sample as the local zero of the experiment.
removeInitialSensorOffsets = true;

% Output weighting during parameter estimation.
outputWeight = diag([1, 8]);

% Estimation settings.
parameterMaxIterations = 140;
icOnlyMaxIterations    = 70;

% Initial-state bounds.
positionMargin = 0.75;    % rad around segment first sample
velocityLimit  = 45.0;    % rad/s

% Actuator-state initial condition bound.
% The actuator state u_act is allowed to start close to the first input,
% but not forced to be exactly equal. This matters when the validation
% segment starts in the middle of a PRBS experiment.
uActMarginFactor = 2.0;
uActMinimumMargin = 0.05;

% Optional: skip the beginning of the training segment if you deliberately
% recorded a pre-input rest period. For your current data, keep this at 0.
skipTrainingStartSeconds = 0.0;

%% ================================================================
% Load and preprocess PRBS data
% ================================================================

raw = load(dataFile);

theta1Name = "theta_1_" + datasetName;
theta2Name = "theta_2_" + datasetName;
uName      = "u_"       + datasetName;

requiredNames = [theta1Name, theta2Name, uName];

for k = 1:numel(requiredNames)
    if ~isfield(raw, requiredNames(k))
        error("Missing variable in data file: %s", requiredNames(k));
    end
end

theta1_ts = raw.(theta1Name);
theta2_ts = raw.(theta2Name);
u_ts      = raw.(uName);

t_raw      = double(theta1_ts.Time(:));
theta1_raw = double(theta1_ts.Data(:));
theta2_raw = double(theta2_ts.Data(:));
u_raw      = double(u_ts.Data(:));

N = min([numel(t_raw), numel(theta1_raw), numel(theta2_raw), numel(u_raw)]);

t_raw      = t_raw(1:N);
theta1_raw = theta1_raw(1:N);
theta2_raw = theta2_raw(1:N);
u_raw      = u_raw(1:N);

t_raw = t_raw - t_raw(1);

idxUsed = t_raw <= usedEndTime;

t      = t_raw(idxUsed);
theta1 = theta1_raw(idxUsed);
phi2   = relativeAngleSign * theta2_raw(idxUsed);
u      = inputPolarity * u_raw(idxUsed);

if removeInitialSensorOffsets
    theta1Offset = theta1(1);
    phi2Offset   = phi2(1);

    theta1 = theta1 - theta1Offset;
    phi2   = phi2   - phi2Offset;
else
    theta1Offset = 0;
    phi2Offset   = 0;
end

t = t - t(1);

Ts = mean(diff(t));
fs = 1 / Ts;

theta2_abs = theta1 + phi2;

fprintf('\nSelected dataset: %s\n', datasetName);
fprintf('Sampling time Ts: %.8f s\n', Ts);
fprintf('Sampling frequency fs: %.3f Hz\n', fs);
fprintf('Used duration: %.3f s\n', t(end));
fprintf('Number of samples used: %d\n', numel(t));
fprintf('Input polarity: %+d\n', inputPolarity);
fprintf('Relative angle sign: %+d\n', relativeAngleSign);
fprintf('Removed theta_1 initial offset: %.8g rad\n', theta1Offset);
fprintf('Removed phi_2 initial offset: %.8g rad\n', phi2Offset);
fprintf('Output weights: theta_1 = %.3g, phi_2 = %.3g\n', outputWeight(1,1), outputWeight(2,2));

plotMeasuredDataAndSpectrum(t, theta1, phi2, theta2_abs, u, fs);

%% ================================================================
% Build training and validation iddata objects
% ================================================================

effectiveTrainStart = trainStart + skipTrainingStartSeconds;

idxTrain = t >= effectiveTrainStart & t < trainEnd;
idxVal   = t >= valStart            & t <= valEnd;

tTrain = t(idxTrain) - t(find(idxTrain, 1, 'first'));
tVal   = t(idxVal)   - t(find(idxVal,   1, 'first'));

yTrain = [theta1(idxTrain), phi2(idxTrain)];
uTrain = u(idxTrain);

yVal = [theta1(idxVal), phi2(idxVal)];
uVal = u(idxVal);

zTrain = iddata(yTrain, uTrain, Ts);
zVal   = iddata(yVal,   uVal,   Ts);

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

fprintf('\nTraining samples:   %d, duration %.3f s\n', numel(uTrain), tTrain(end));
fprintf('Validation samples: %d, duration %.3f s\n', numel(uVal), tVal(end));

% Input margins for actuator-state IC.
uActMarginTrain = max(uActMinimumMargin, uActMarginFactor * max(abs(uTrain)));
uActMarginVal   = max(uActMinimumMargin, uActMarginFactor * max(abs(uVal)));

%% ================================================================
% Create upgraded 8-parameter white-box model
% ================================================================

modelFile = 'final_whitebox_8par_lag_bias_model';

Order = [2 1 5];  % 2 outputs, 1 input, 5 states
TsModel = 0;      % continuous-time model

% Initial guesses based on previous successful 6-parameter run.
p0.p_J1    = 14.0;
p0.p_kappa = 0.31;
p0.p_g1    = 97.0;
p0.p_g2    = 112.0;
p0.p_Ku    = 4020.0;
p0.p_b1    = 670.0;

% New parameters:
%   p_tau0 models constant torque/input bias.
%   p_Tm models first-order actuator lag.
p0.p_tau0  = 0.0;
p0.p_Tm    = 0.02;

theta1_0_train = yTrain(1, 1);
phi2_0_train   = yTrain(1, 2);
theta2_0_train = theta1_0_train + phi2_0_train;
uact_0_train   = uTrain(1);

InitialStates = {theta1_0_train; theta2_0_train; 0; 0; uact_0_train};

Parameters = {
    p0.p_J1
    p0.p_kappa
    p0.p_g1
    p0.p_g2
    p0.p_Ku
    p0.p_b1
    p0.p_tau0
    p0.p_Tm
};

model0 = idnlgrey(modelFile, Order, Parameters, InitialStates, TsModel);
model0 = configureFinal8ParModel(model0);

% Estimate initial conditions during parameter identification.
model0 = configureInitialStateBounds(model0, yTrain, uTrain, ...
    positionMargin, velocityLimit, uActMarginTrain, 'training');

fprintf('\nInitial model:\n');
dispParameters(model0);
dispInitialStates(model0);

%% ================================================================
% Estimate physical parameters using training data
% ================================================================

optParam = nlgreyestOptions;
optParam.Display = 'on';
optParam.EstimateCovariance = true;
optParam.SearchOptions.MaxIterations = parameterMaxIterations;
optParam.OutputWeight = outputWeight;

fprintf('\n============================================================\n');
fprintf('STEP 1: Estimate 8 parameters and training initial states\n');
fprintf('============================================================\n');

modelParam = nlgreyest(zTrain, model0, optParam);
modelParam.Name = 'Final 8-parameter lag/bias white-box model';

fprintf('\nEstimated parameters after training fit:\n');
dispParameters(modelParam);

fprintf('Estimated training initial states after parameter fit:\n');
dispInitialStates(modelParam);

%% ================================================================
% Re-estimate IC only for training and validation, parameters fixed
% ================================================================

optIC = nlgreyestOptions;
optIC.Display = 'on';
optIC.EstimateCovariance = false;
optIC.SearchOptions.MaxIterations = icOnlyMaxIterations;
optIC.OutputWeight = outputWeight;

fprintf('\n============================================================\n');
fprintf('STEP 2: Estimate training IC only, parameters fixed\n');
fprintf('============================================================\n');

[modelTrainIC, fitTrain] = estimateICOnlyAndCompare( ...
    zTrain, modelParam, yTrain, uTrain, positionMargin, velocityLimit, ...
    uActMarginTrain, 'training', optIC);

fprintf('\nTraining IC-only result:\n');
dispInitialStates(modelTrainIC);
dispFitBreakdown(fitTrain, {'theta_1', 'phi_2'});

fprintf('\n============================================================\n');
fprintf('STEP 3: Estimate validation IC only, parameters fixed\n');
fprintf('============================================================\n');

[modelValIC, fitVal] = estimateICOnlyAndCompare( ...
    zVal, modelParam, yVal, uVal, positionMargin, velocityLimit, ...
    uActMarginVal, 'validation', optIC);

fprintf('\nValidation IC-only result:\n');
dispInitialStates(modelValIC);
dispFitBreakdown(fitVal, {'theta_1', 'phi_2'});

%% ================================================================
% Final plots and diagnostics
% ================================================================

[yHatTrain, fitTrain] = compareUsingModelIC(zTrain, modelTrainIC);
[yHatVal,   fitVal]   = compareUsingModelIC(zVal,   modelValIC);

fprintf('\n============================================================\n');
fprintf('FINAL FIT WITH SEGMENT-SPECIFIC ESTIMATED ICs\n');
fprintf('============================================================\n');

fprintf('\nTraining fit:\n');
dispFitBreakdown(fitTrain, {'theta_1', 'phi_2'});

fprintf('Validation fit:\n');
dispFitBreakdown(fitVal, {'theta_1', 'phi_2'});

fprintf('\nFinal parameter values:\n');
dispParameters(modelParam);

fprintf('\nTraining initial states:\n');
dispInitialStates(modelTrainIC);

fprintf('Validation initial states:\n');
dispInitialStates(modelValIC);

plotFinalComparison(zTrain, zVal, yHatTrain, yHatVal, tTrain, tVal, fitTrain, fitVal);
plotResidualDiagnostics(zTrain, zVal, yHatTrain, yHatVal, tTrain, tVal);

fprintf('\nAmplitude diagnostics, training:\n');
printAmplitudeDiagnostics(zTrain, yHatTrain, {'theta_1', 'phi_2'});

fprintf('\nAmplitude diagnostics, validation:\n');
printAmplitudeDiagnostics(zVal, yHatVal, {'theta_1', 'phi_2'});

fprintf('\nResidual diagnostics, training:\n');
printResidualDiagnostics(zTrain, yHatTrain, {'theta_1', 'phi_2'});

fprintf('\nResidual diagnostics, validation:\n');
printResidualDiagnostics(zVal, yHatVal, {'theta_1', 'phi_2'});

fprintf('\nInterpretation of new parameters:\n');
fprintf('  p_tau0 = %.8g. If this is far from 0, the previous offset was partly a torque/input bias.\n', ...
    getParValue(modelParam, 'p_tau0'));
fprintf('  p_Tm   = %.8g s. If this is near its lower bound, actuator lag is not important.\n', ...
    getParValue(modelParam, 'p_Tm'));

%% ================================================================
% Save result
% ================================================================

finalWhitebox8.modelParam = modelParam;
finalWhitebox8.modelTrainIC = modelTrainIC;
finalWhitebox8.modelValIC = modelValIC;
finalWhitebox8.fitTrain = fitTrain;
finalWhitebox8.fitVal = fitVal;

finalWhitebox8.settings.dataFile = dataFile;
finalWhitebox8.settings.datasetName = datasetName;
finalWhitebox8.settings.inputPolarity = inputPolarity;
finalWhitebox8.settings.relativeAngleSign = relativeAngleSign;
finalWhitebox8.settings.removeInitialSensorOffsets = removeInitialSensorOffsets;
finalWhitebox8.settings.theta1Offset = theta1Offset;
finalWhitebox8.settings.phi2Offset = phi2Offset;
finalWhitebox8.settings.trainStart = trainStart;
finalWhitebox8.settings.trainEnd = trainEnd;
finalWhitebox8.settings.valStart = valStart;
finalWhitebox8.settings.valEnd = valEnd;
finalWhitebox8.settings.skipTrainingStartSeconds = skipTrainingStartSeconds;
finalWhitebox8.settings.outputWeight = outputWeight;
finalWhitebox8.settings.positionMargin = positionMargin;
finalWhitebox8.settings.velocityLimit = velocityLimit;
finalWhitebox8.settings.uActMarginTrain = uActMarginTrain;
finalWhitebox8.settings.uActMarginVal = uActMarginVal;

save('final_whitebox_8par_lag_bias_identification_result.mat', 'finalWhitebox8');

fprintf('\nSaved result to final_whitebox_8par_lag_bias_identification_result.mat\n');

%% ================================================================
% Local helper functions
% ================================================================

function model = configureFinal8ParModel(model)

names = {
    'p_J1'
    'p_kappa'
    'p_g1'
    'p_g2'
    'p_Ku'
    'p_b1'
    'p_tau0'
    'p_Tm'
};

units = {
    '-'
    '-'
    '1/s^2'
    '1/s^2'
    'torque/input'
    '1/s'
    'torque'
    's'
};

mins = [
    0.05
   -0.98
    0
    0
    0
    0
   -3000
    0.001
];

maxs = [
    60
    0.98
    700
    700
    15000
    4000
    3000
    0.50
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

% MATLAB R2021b idnlgrey state metadata.
model.InitialStates(1).Name = 'theta_1';
model.InitialStates(1).Unit = 'rad';
model.InitialStates(2).Name = 'theta_2_abs';
model.InitialStates(2).Unit = 'rad';
model.InitialStates(3).Name = 'omega_1';
model.InitialStates(3).Unit = 'rad/s';
model.InitialStates(4).Name = 'omega_2_abs';
model.InitialStates(4).Unit = 'rad/s';
model.InitialStates(5).Name = 'u_act';
model.InitialStates(5).Unit = 'input';

end

function model = configureInitialStateBounds(model, y, uSegment, ...
    positionMargin, velocityLimit, uActMargin, label)

theta1_0 = y(1, 1);
phi2_0   = y(1, 2);
theta2_0 = theta1_0 + phi2_0;

uact_0 = uSegment(1);

values = {theta1_0; theta2_0; 0; 0; uact_0};
mins   = {theta1_0 - positionMargin; theta2_0 - positionMargin; -velocityLimit; -velocityLimit; uact_0 - uActMargin};
maxs   = {theta1_0 + positionMargin; theta2_0 + positionMargin;  velocityLimit;  velocityLimit; uact_0 + uActMargin};
fixed  = {false; false; false; false; false};

% Reset to wide bounds first. MATLAB checks whether the current value lies
% inside the bounds whenever setinit is called, so this avoids errors when
% switching from training ICs to validation ICs.
wideMin = {-Inf; -Inf; -Inf; -Inf; -Inf};
wideMax = { Inf;  Inf;  Inf;  Inf;  Inf};

model = setinit(model, 'Minimum', wideMin);
model = setinit(model, 'Maximum', wideMax);
model = setinit(model, 'Value', values);
model = setinit(model, 'Minimum', mins);
model = setinit(model, 'Maximum', maxs);
model = setinit(model, 'Fixed', fixed);

fprintf('Configured %s IC bounds around first measured sample and first input.\n', label);

end

function [modelIC, fit] = estimateICOnlyAndCompare(z, modelParam, y, uSegment, ...
    positionMargin, velocityLimit, uActMargin, label, optIC)

modelIC = modelParam;

% Fix all physical parameters. Only the initial states may move.
for k = 1:length(modelIC.Parameters)
    modelIC.Parameters(k).Fixed = true;
end

modelIC = configureInitialStateBounds(modelIC, y, uSegment, ...
    positionMargin, velocityLimit, uActMargin, label);

modelIC = nlgreyest(z, modelIC, optIC);
modelIC.Name = sprintf('%s IC-only model', label);

[~, fit] = compareUsingModelIC(z, modelIC);

end

function [yHat, fit] = compareUsingModelIC(data, model)

% Prefer comparing with the model's own stored initial condition.
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

function plotMeasuredDataAndSpectrum(t, theta1, phi2, theta2_abs, u, fs)

figure('Name', 'Measured PRBS data, first 90 seconds', ...
    'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.8]);
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t, theta1);
grid on;
ylabel('\theta_1 [rad]');
title('Measured rod-1 absolute angle');

nexttile;
plot(t, phi2);
grid on;
ylabel('\phi_2 [rad]');
title('Measured rod-2 relative encoder angle');

nexttile;
plot(t, theta2_abs);
grid on;
ylabel('\theta_{2,abs} [rad]');
title('Computed absolute rod-2 angle: \theta_{2,abs} = \theta_1 + \phi_2');

nexttile;
plot(t, u);
grid on;
ylabel('u');
xlabel('Time [s]');
title('Input');

figure('Name', 'Fourier spectra', ...
    'Units', 'normalized', 'Position', [0.08 0.08 0.85 0.75]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

plotSpectrumTile(theta1, fs, '\theta_1 spectrum');
plotSpectrumTile(phi2,   fs, '\phi_2 spectrum');
plotSpectrumTile(u,      fs, 'input spectrum');

end

function plotSpectrumTile(signal, fs, plotTitle)

nexttile;

signal = signal(:);
signal = signal - mean(signal, 'omitnan');

N = numel(signal);
Y = fft(signal);
P2 = abs(Y / N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);

f = fs * (0:floor(N/2)) / N;

plot(f, P1);
grid on;
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title(plotTitle);
xlim([0, min(20, fs/2)]);

end

function plotFinalComparison(zTrain, zVal, yHatTrain, yHatVal, tTrain, tVal, fitTrain, fitVal)

yTrain = getOutputDataMatrix(zTrain);
yVal   = getOutputDataMatrix(zVal);

yTrainHat = getOutputDataMatrix(yHatTrain);
yValHat   = getOutputDataMatrix(yHatVal);

trainFitMean = meanFitValue(fitTrain);
valFitMean   = meanFitValue(fitVal);

figure('Name', 'Final 8-parameter white-box lag/bias fit', ...
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

sgtitle('Final 8-parameter white-box model with actuator lag/bias and segment-specific ICs');

end

function plotResidualDiagnostics(zTrain, zVal, yHatTrain, yHatVal, tTrain, tVal)

yTrain = getOutputDataMatrix(zTrain);
yVal   = getOutputDataMatrix(zVal);

yTrainHat = getOutputDataMatrix(yHatTrain);
yValHat   = getOutputDataMatrix(yHatVal);

nTrain = min(size(yTrain,1), size(yTrainHat,1));
nVal   = min(size(yVal,1),   size(yValHat,1));

rTrain = yTrain(1:nTrain,:) - yTrainHat(1:nTrain,:);
rVal   = yVal(1:nVal,:)     - yValHat(1:nVal,:);

figure('Name', 'Residual diagnostics', ...
    'Units', 'normalized', 'Position', [0.08 0.08 0.85 0.75]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(tTrain(1:nTrain), rTrain(:,1));
grid on;
ylabel('\theta_1 residual [rad]');
title('Training residual \theta_1');

nexttile;
plot(tTrain(1:nTrain), rTrain(:,2));
grid on;
ylabel('\phi_2 residual [rad]');
title('Training residual \phi_2');

nexttile;
plot(tVal(1:nVal), rVal(:,1));
grid on;
ylabel('\theta_1 residual [rad]');
xlabel('Time [s]');
title('Validation residual \theta_1');

nexttile;
plot(tVal(1:nVal), rVal(:,2));
grid on;
ylabel('\phi_2 residual [rad]');
xlabel('Time [s]');
title('Validation residual \phi_2');

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

function printResidualDiagnostics(dataMeasured, dataModel, outputNames)

ym = getOutputDataMatrix(dataMeasured);
yh = getOutputDataMatrix(dataModel);

n = min(size(ym,1), size(yh,1));
ym = ym(1:n,:);
yh = yh(1:n,:);

r = ym - yh;

nOut = min(size(r,2), numel(outputNames));

for k = 1:nOut
    fprintf('  %-8s mean = %.6g, std = %.6g, RMS = %.6g\n', ...
        outputNames{k}, mean(r(:,k)), std(r(:,k)), sqrt(mean(r(:,k).^2)));
end

end

function value = getParValue(model, parName)

for k = 1:length(model.Parameters)
    if strcmp(model.Parameters(k).Name, parName)
        value = double(model.Parameters(k).Value);
        return;
    end
end

error('Parameter %s not found.', parName);

end
