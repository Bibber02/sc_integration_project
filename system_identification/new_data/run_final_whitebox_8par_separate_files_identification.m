clear;
clear functions;
clc;
close all;

fprintf('Running run_final_whitebox_8par_separate_files_identification.m\n');

%% ================================================================
% User settings
% ================================================================

% New separate recordings.
% The script first tries these paths exactly as written. If it cannot find
% the files, it also tries the same names inside dataFolder.
dataFolder = "identification_data\whole_system_identification";
trainingFileName   = "training_data_amp_0125_120s.mat";
validationFileName = "validation_data_amp_0_125_120s.mat";

trainingFile   = resolveDataFile(trainingFileName, dataFolder);
validationFile = resolveDataFile(validationFileName, dataFolder);

% Variables inside both files:
%   u_ts, theta_1, theta_2
inputVariableName  = "u_ts";
theta1VariableName = "theta_1";
theta2VariableName = "theta_2";

% Use the full recordings by default.
% Set to e.g. 2.0 if you recorded a zero-input settling period and want to
% exclude it from parameter estimation. Usually keep it at 0 first.
trainingStartTime   = 0.0;
validationStartTime = 0.0;

% Use Inf to take the full file duration.
trainingEndTime   = Inf;
validationEndTime = Inf;

% Sign conventions from earlier tests.
inputPolarity = -1;
relativeAngleSign = 1;

% Treat the first sample of each recording as the local zero of that recording.
% This is usually convenient when encoders have arbitrary zero offsets.
removeInitialSensorOffsets = true;

% Output weighting during parameter estimation.
% Leave phi_2 high because theta_1 is usually easier to fit.
outputWeight = diag([1, 8]);

% Estimation settings.
parameterMaxIterations = 150;
icOnlyMaxIterations    = 70;

% Initial-state bounds.
positionMargin = 0.75;    % rad around segment first sample
velocityLimit  = 45.0;    % rad/s

% Actuator-state IC bounds. If your input is normalized PRBS around +/-0.125,
% these values are intentionally broad.
uActMargin = 1.0;

% Optional decimation for faster early testing.
% Set to 1 for final identification. Set to 2 or 4 to test quickly.
decimationFactor = 1;

%% ================================================================
% Load and preprocess training and validation recordings
% ================================================================

trainRaw = loadRecording(trainingFile, inputVariableName, theta1VariableName, theta2VariableName);
valRaw   = loadRecording(validationFile, inputVariableName, theta1VariableName, theta2VariableName);

trainData = preprocessRecording(trainRaw, inputPolarity, relativeAngleSign, ...
    removeInitialSensorOffsets, trainingStartTime, trainingEndTime, decimationFactor, 'training');

valData = preprocessRecording(valRaw, inputPolarity, relativeAngleSign, ...
    removeInitialSensorOffsets, validationStartTime, validationEndTime, decimationFactor, 'validation');

fprintf('\nTraining file:   %s\n', trainingFile);
fprintf('Validation file: %s\n', validationFile);
fprintf('Input polarity: %+d\n', inputPolarity);
fprintf('Relative angle sign: %+d\n', relativeAngleSign);
fprintf('Output weights: theta_1 = %.3g, phi_2 = %.3g\n', outputWeight(1,1), outputWeight(2,2));
fprintf('Decimation factor: %d\n', decimationFactor);

fprintf('\nTraining data:\n');
printDataSummary(trainData);

fprintf('\nValidation data:\n');
printDataSummary(valData);

plotRecordingOverview(trainData, valData);

%% ================================================================
% Build iddata objects
% ================================================================

zTrain = iddata([trainData.theta1, trainData.phi2], trainData.u, trainData.Ts);
zVal   = iddata([valData.theta1,   valData.phi2],   valData.u,   valData.Ts);

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
% Create 8-parameter white-box model with actuator lag and bias
% ================================================================

modelFile = 'final_whitebox_8par_lag_bias_model_separate_files';

Order = [2 1 5];  % 2 outputs, 1 input, 5 states
TsModel = 0;      % continuous-time model

% Initial guesses based on earlier successful runs.
p0.p_J1    = 14.0;
p0.p_kappa = 0.30;
p0.p_g1    = 97.0;
p0.p_g2    = 112.0;
p0.p_Ku    = 4000.0;
p0.p_b1    = 670.0;
p0.p_tau0  = 0.0;
p0.p_Tm    = 0.02;

InitialStates = makeInitialStatesFromData(trainData);

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
model0 = configure8ParModel(model0);

% Estimate initial conditions during parameter identification.
model0 = configureInitialStateBounds(model0, trainData, positionMargin, velocityLimit, uActMargin, 'training');

fprintf('\nInitial model:\n');
dispParameters(model0);
dispInitialStates(model0);

%% ================================================================
% Estimate physical parameters using training recording
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
modelParam.Name = 'Final 8-parameter white-box model, separate files';

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
    zTrain, modelParam, trainData, positionMargin, velocityLimit, uActMargin, ...
    'training', optIC);

fprintf('\nTraining IC-only result:\n');
dispInitialStates(modelTrainIC);
dispFitBreakdown(fitTrain, {'theta_1', 'phi_2'});

fprintf('\n============================================================\n');
fprintf('STEP 3: Estimate validation IC only, parameters fixed\n');
fprintf('============================================================\n');

[modelValIC, fitVal] = estimateICOnlyAndCompare( ...
    zVal, modelParam, valData, positionMargin, velocityLimit, uActMargin, ...
    'validation', optIC);

fprintf('\nValidation IC-only result:\n');
dispInitialStates(modelValIC);
dispFitBreakdown(fitVal, {'theta_1', 'phi_2'});

%% ================================================================
% Final plots and diagnostics
% ================================================================

[yHatTrain, fitTrain] = compareUsingModelIC(zTrain, modelTrainIC);
[yHatVal,   fitVal]   = compareUsingModelIC(zVal,   modelValIC);

fprintf('\n============================================================\n');
fprintf('FINAL FIT WITH SEPARATE TRAINING/VALIDATION FILES\n');
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

plotFinalComparison(zTrain, zVal, yHatTrain, yHatVal, trainData.t, valData.t, fitTrain, fitVal);

fprintf('\nAmplitude diagnostics, training:\n');
printAmplitudeDiagnostics(zTrain, yHatTrain, {'theta_1', 'phi_2'});

fprintf('\nAmplitude diagnostics, validation:\n');
printAmplitudeDiagnostics(zVal, yHatVal, {'theta_1', 'phi_2'});

%% ================================================================
% Save result
% ================================================================

finalWhiteboxSeparate.modelParam = modelParam;
finalWhiteboxSeparate.modelTrainIC = modelTrainIC;
finalWhiteboxSeparate.modelValIC = modelValIC;
finalWhiteboxSeparate.fitTrain = fitTrain;
finalWhiteboxSeparate.fitVal = fitVal;
finalWhiteboxSeparate.trainingFile = trainingFile;
finalWhiteboxSeparate.validationFile = validationFile;
finalWhiteboxSeparate.trainData = trainData;
finalWhiteboxSeparate.valData = valData;

finalWhiteboxSeparate.settings.inputPolarity = inputPolarity;
finalWhiteboxSeparate.settings.relativeAngleSign = relativeAngleSign;
finalWhiteboxSeparate.settings.removeInitialSensorOffsets = removeInitialSensorOffsets;
finalWhiteboxSeparate.settings.trainingStartTime = trainingStartTime;
finalWhiteboxSeparate.settings.validationStartTime = validationStartTime;
finalWhiteboxSeparate.settings.trainingEndTime = trainingEndTime;
finalWhiteboxSeparate.settings.validationEndTime = validationEndTime;
finalWhiteboxSeparate.settings.outputWeight = outputWeight;
finalWhiteboxSeparate.settings.positionMargin = positionMargin;
finalWhiteboxSeparate.settings.velocityLimit = velocityLimit;
finalWhiteboxSeparate.settings.uActMargin = uActMargin;
finalWhiteboxSeparate.settings.decimationFactor = decimationFactor;

save('final_whitebox_8par_separate_files_result.mat', 'finalWhiteboxSeparate');

fprintf('\nSaved result to final_whitebox_8par_separate_files_result.mat\n');

%% ================================================================
% Local helper functions
% ================================================================

function fullPath = resolveDataFile(fileName, dataFolder)

if isfile(fileName)
    fullPath = fileName;
    return;
end

candidate = fullfile(dataFolder, fileName);

if isfile(candidate)
    fullPath = candidate;
    return;
end

error('Could not find data file "%s" either in the current folder or in "%s".', fileName, dataFolder);

end

function rec = loadRecording(fileName, inputName, theta1Name, theta2Name)

S = load(fileName);

if ~isfield(S, inputName)
    error('File %s does not contain variable %s.', fileName, inputName);
end

if ~isfield(S, theta1Name)
    error('File %s does not contain variable %s.', fileName, theta1Name);
end

if ~isfield(S, theta2Name)
    error('File %s does not contain variable %s.', fileName, theta2Name);
end

[tU, u]       = extractSignal(S.(inputName), inputName);
[t1, theta1]  = extractSignal(S.(theta1Name), theta1Name);
[t2, theta2]  = extractSignal(S.(theta2Name), theta2Name);

% Choose a common time base. Prefer theta_1 time, then input time, then theta_2.
if ~isempty(t1)
    t = t1;
elseif ~isempty(tU)
    t = tU;
elseif ~isempty(t2)
    t = t2;
else
    % Fallback only if all signals are plain vectors without time.
    % Your earlier data used 100 Hz, so Ts = 0.01 is the safest default.
    TsFallback = 0.01;
    N = min([numel(u), numel(theta1), numel(theta2)]);
    t = (0:N-1)' * TsFallback;
end

t = double(t(:));
t = t - t(1);

u      = resampleIfNeeded(tU, u, t);
theta1 = resampleIfNeeded(t1, theta1, t);
theta2 = resampleIfNeeded(t2, theta2, t);

N = min([numel(t), numel(u), numel(theta1), numel(theta2)]);

rec.fileName = fileName;
rec.t = t(1:N);
rec.u = u(1:N);
rec.theta1 = theta1(1:N);
rec.theta2 = theta2(1:N);

end

function [t, data] = extractSignal(obj, name)

if isa(obj, 'timeseries')
    t = double(obj.Time(:));
    data = double(obj.Data(:));
    return;
end

if isstruct(obj) && isfield(obj, 'Time') && isfield(obj, 'Data')
    t = double(obj.Time(:));
    data = double(obj.Data(:));
    return;
end

% Some Simulink logged signals are stored as a structure with signals.values.
if isstruct(obj) && isfield(obj, 'time') && isfield(obj, 'signals')
    t = double(obj.time(:));
    data = double(obj.signals.values(:));
    return;
end

if isnumeric(obj)
    t = [];
    data = double(obj(:));
    return;
end

error('Unsupported format for variable %s. Expected timeseries, struct with Time/Data, or numeric vector.', name);

end

function y = resampleIfNeeded(tOriginal, yOriginal, tTarget)

yOriginal = double(yOriginal(:));

if isempty(tOriginal)
    % No time vector. Assume it is already sampled on target time.
    N = min(numel(yOriginal), numel(tTarget));
    y = yOriginal(1:N);

    if N < numel(tTarget)
        y(end+1:numel(tTarget),1) = y(end);
    end

    return;
end

tOriginal = double(tOriginal(:));
tOriginal = tOriginal - tOriginal(1);

if numel(tOriginal) == numel(tTarget) && max(abs(tOriginal - tTarget)) < 1e-10
    y = yOriginal;
else
    y = interp1(tOriginal, yOriginal, tTarget, 'linear', 'extrap');
end

y = double(y(:));

end

function data = preprocessRecording(raw, inputPolarity, relativeAngleSign, removeOffsets, tStart, tEnd, decimationFactor, label)

idx = raw.t >= tStart;

if isfinite(tEnd)
    idx = idx & raw.t <= tEnd;
end

t = raw.t(idx);
theta1 = raw.theta1(idx);
phi2 = relativeAngleSign * raw.theta2(idx);
u = inputPolarity * raw.u(idx);

% Shift selected time window to start at zero.
t = t - t(1);

if removeOffsets
    theta1Offset = theta1(1);
    phi2Offset   = phi2(1);

    theta1 = theta1 - theta1Offset;
    phi2   = phi2   - phi2Offset;
else
    theta1Offset = 0;
    phi2Offset   = 0;
end

if decimationFactor > 1
    t = t(1:decimationFactor:end);
    theta1 = theta1(1:decimationFactor:end);
    phi2 = phi2(1:decimationFactor:end);
    u = u(1:decimationFactor:end);
end

Ts = mean(diff(t));

data.label = label;
data.fileName = raw.fileName;
data.t = t(:);
data.theta1 = theta1(:);
data.phi2 = phi2(:);
data.theta2_abs = data.theta1 + data.phi2;
data.u = u(:);
data.Ts = Ts;
data.fs = 1/Ts;
data.theta1Offset = theta1Offset;
data.phi2Offset = phi2Offset;

end

function printDataSummary(data)

fprintf('  File: %s\n', data.fileName);
fprintf('  Samples: %d\n', numel(data.t));
fprintf('  Duration: %.3f s\n', data.t(end));
fprintf('  Ts: %.8f s\n', data.Ts);
fprintf('  fs: %.3f Hz\n', data.fs);
fprintf('  Removed theta_1 offset: %.8g rad\n', data.theta1Offset);
fprintf('  Removed phi_2 offset: %.8g rad\n', data.phi2Offset);
fprintf('  Input min/max: %.6g / %.6g\n', min(data.u), max(data.u));

end

function InitialStates = makeInitialStatesFromData(data)

theta1_0 = data.theta1(1);
phi2_0   = data.phi2(1);
theta2_0 = theta1_0 + phi2_0;
uAct_0 = data.u(1);

InitialStates = {theta1_0; theta2_0; 0; 0; uAct_0};

end

function model = configure8ParModel(model)

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
   -5000
    0.005
];

maxs = [
    50
    0.98
    600
    600
    15000
    4000
    5000
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

function model = configureInitialStateBounds(model, data, positionMargin, velocityLimit, uActMargin, label)

theta1_0 = data.theta1(1);
phi2_0   = data.phi2(1);
theta2_0 = theta1_0 + phi2_0;
uAct_0 = data.u(1);

values = {theta1_0; theta2_0; 0; 0; uAct_0};
mins   = {theta1_0 - positionMargin; theta2_0 - positionMargin; -velocityLimit; -velocityLimit; uAct_0 - uActMargin};
maxs   = {theta1_0 + positionMargin; theta2_0 + positionMargin;  velocityLimit;  velocityLimit; uAct_0 + uActMargin};
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

fprintf('Configured %s IC bounds around first measured sample.\n', label);

end

function [modelIC, fit] = estimateICOnlyAndCompare(z, modelParam, data, ...
    positionMargin, velocityLimit, uActMargin, label, optIC)

modelIC = modelParam;

% Fix all physical parameters. Only the initial states may move.
for k = 1:length(modelIC.Parameters)
    modelIC.Parameters(k).Fixed = true;
end

modelIC = configureInitialStateBounds(modelIC, data, positionMargin, velocityLimit, uActMargin, label);

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

function plotRecordingOverview(trainData, valData)

figure('Name', 'Measured separate training and validation recordings', ...
    'Units', 'normalized', 'Position', [0.04 0.05 0.9 0.85]);
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(trainData.t, trainData.theta1);
grid on;
ylabel('\theta_1 [rad]');
title('Training \theta_1');

nexttile;
plot(valData.t, valData.theta1);
grid on;
title('Validation \theta_1');

nexttile;
plot(trainData.t, trainData.phi2);
grid on;
ylabel('\phi_2 [rad]');
title('Training \phi_2');

nexttile;
plot(valData.t, valData.phi2);
grid on;
title('Validation \phi_2');

nexttile;
plot(trainData.t, trainData.theta2_abs);
grid on;
ylabel('\theta_{2,abs} [rad]');
title('Training computed \theta_{2,abs}');

nexttile;
plot(valData.t, valData.theta2_abs);
grid on;
title('Validation computed \theta_{2,abs}');

nexttile;
plot(trainData.t, trainData.u);
grid on;
ylabel('u');
xlabel('Time [s]');
title('Training input');

nexttile;
plot(valData.t, valData.u);
grid on;
xlabel('Time [s]');
title('Validation input');

figure('Name', 'Input and output spectra', ...
    'Units', 'normalized', 'Position', [0.08 0.08 0.85 0.75]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

plotSpectrumTile(trainData.theta1, trainData.fs, 'Training \theta_1 spectrum');
plotSpectrumTile(valData.theta1,   valData.fs,   'Validation \theta_1 spectrum');
plotSpectrumTile(trainData.phi2,   trainData.fs, 'Training \phi_2 spectrum');
plotSpectrumTile(valData.phi2,     valData.fs,   'Validation \phi_2 spectrum');
plotSpectrumTile(trainData.u,      trainData.fs, 'Training input spectrum');
plotSpectrumTile(valData.u,        valData.fs,   'Validation input spectrum');

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

figure('Name', 'Final 8-parameter white-box fit, separate files', ...
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

sgtitle('Final 8-parameter white-box model, separate training and validation recordings');

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
