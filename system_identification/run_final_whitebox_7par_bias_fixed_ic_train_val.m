clear;
clear functions;
clc;
close all;

fprintf('Running run_final_whitebox_7par_bias_fixed_ic_train_val.m\n');

%% ================================================================
% User settings
% ================================================================

% New file with two separate PRBS recordings.
% Change this to the actual location/name of your new .mat file.
dataFile = "identification_data\whole_system_identification\prbs_train_validation_121s.mat";

% The script expects variables named:
%   theta_1_<datasetName>
%   theta_2_<datasetName>
%   u_<datasetName>
%
% Example:
%   trainDatasetName      = "Train";
%   validationDatasetName = "Validation";
% expects:
%   theta_1_Train, theta_2_Train, u_Train
%   theta_1_Validation, theta_2_Validation, u_Validation
%
% If your variable names are different, fill in the direct names below.
trainDatasetName      = "Train";
validationDatasetName = "Validation";

% Direct variable-name override. Leave these as "" to use the pattern above.
trainVar.theta1 = "";
trainVar.theta2 = "";
trainVar.u      = "";

valVar.theta1 = "";
valVar.theta2 = "";
valVar.u      = "";

% You recorded 121 seconds with the first second idle.
usedStartTime = 0.0;
usedEndTime   = 121.0;

% Keep the idle second in the estimation by default.
% If you later want to estimate only on the PRBS portion, set this to 1.0.
estimationStartTime = 0.0;

% Sign conventions from previous successful identification.
inputPolarity = -1;
relativeAngleSign = 1;

% If the first second is a rest/zero reference, this is usually correct.
% It sets the initial measured angles to zero for each run separately.
removeInitialSensorOffsets = true;

% Output weighting.
outputWeight = diag([1, 8]);

% Estimation settings.
parameterMaxIterations = 150;

% Optional second parameter-only refinement.
runSecondRefinement = true;
secondRefinementIterations = 80;

% Optional multistart. Keep false unless the first result is poor.
runMultistart = false;
nMultistarts = 4;
multistartIterations = 80;

%% ================================================================
% Load and preprocess data
% ================================================================

raw = load(dataFile);

[trainData, trainInfo] = loadAndPreprocessRun( ...
    raw, trainDatasetName, trainVar, usedStartTime, usedEndTime, ...
    inputPolarity, relativeAngleSign, removeInitialSensorOffsets, "training");

[valData, valInfo] = loadAndPreprocessRun( ...
    raw, validationDatasetName, valVar, usedStartTime, usedEndTime, ...
    inputPolarity, relativeAngleSign, removeInitialSensorOffsets, "validation");

TsTrain = mean(diff(trainData.t));
TsVal   = mean(diff(valData.t));

if abs(TsTrain - TsVal) > 1e-8
    warning('Training and validation sampling times differ: %.12g vs %.12g', TsTrain, TsVal);
end

Ts = TsTrain;
fs = 1 / Ts;

fprintf('\nLoaded training dataset: %s\n', trainInfo.name);
fprintf('  variables: %s, %s, %s\n', trainInfo.theta1Name, trainInfo.theta2Name, trainInfo.uName);
fprintf('  duration: %.3f s, samples: %d, Ts: %.8f s\n', trainData.t(end), numel(trainData.t), TsTrain);
fprintf('  theta_1 offset removed: %.8g rad\n', trainInfo.theta1Offset);
fprintf('  phi_2 offset removed:   %.8g rad\n', trainInfo.phi2Offset);

fprintf('\nLoaded validation dataset: %s\n', valInfo.name);
fprintf('  variables: %s, %s, %s\n', valInfo.theta1Name, valInfo.theta2Name, valInfo.uName);
fprintf('  duration: %.3f s, samples: %d, Ts: %.8f s\n', valData.t(end), numel(valData.t), TsVal);
fprintf('  theta_1 offset removed: %.8g rad\n', valInfo.theta1Offset);
fprintf('  phi_2 offset removed:   %.8g rad\n', valInfo.phi2Offset);

fprintf('\nInput polarity: %+d\n', inputPolarity);
fprintf('Relative angle sign: %+d\n', relativeAngleSign);
fprintf('Output weights: theta_1 = %.3g, phi_2 = %.3g\n', outputWeight(1,1), outputWeight(2,2));

plotMeasuredRuns(trainData, valData, fs);

%% ================================================================
% Build iddata objects
% ================================================================

idxTrain = trainData.t >= estimationStartTime;
idxVal   = valData.t   >= estimationStartTime;

tTrain = trainData.t(idxTrain) - trainData.t(find(idxTrain, 1, 'first'));
tVal   = valData.t(idxVal)     - valData.t(find(idxVal,   1, 'first'));

yTrain = [trainData.theta1(idxTrain), trainData.phi2(idxTrain)];
uTrain = trainData.u(idxTrain);

yVal = [valData.theta1(idxVal), valData.phi2(idxVal)];
uVal = valData.u(idxVal);

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

fprintf('\nEstimation data:\n');
fprintf('  training samples:   %d, duration %.3f s\n', numel(uTrain), tTrain(end));
fprintf('  validation samples: %d, duration %.3f s\n', numel(uVal), tVal(end));

%% ================================================================
% Create model with fixed initial states
% ================================================================

modelFile = 'final_whitebox_7par_bias_fixed_ic_model';
Order = [2 1 4];  % 2 outputs, 1 input, 4 states
TsModel = 0;      % continuous-time model

% Initial guesses based on the previous successful result.
p0.p_J1    = 14.6;
p0.p_kappa = 0.32;
p0.p_g1    = 96.0;
p0.p_g2    = 112.0;
p0.p_Ku    = 4243.0;
p0.p_b1    = 755.0;
p0.p_tau0  = -4.2;

% Fixed initial state from the first sample of the estimation segment.
% With the first second idle and offsets removed, this should be close to zero.
x0Train = fixedInitialStateFromFirstSample(yTrain);
x0Val   = fixedInitialStateFromFirstSample(yVal);

Parameters = {
    p0.p_J1
    p0.p_kappa
    p0.p_g1
    p0.p_g2
    p0.p_Ku
    p0.p_b1
    p0.p_tau0
};

model0 = idnlgrey(modelFile, Order, Parameters, x0Train, TsModel);
model0 = configure7ParModel(model0);

fprintf('\nInitial model:\n');
dispParameters(model0);
dispInitialStates(model0);

%% ================================================================
% Estimate parameters on training data, fixed IC
% ================================================================

optParam = nlgreyestOptions;
optParam.Display = 'on';
optParam.EstimateCovariance = true;
optParam.SearchOptions.MaxIterations = parameterMaxIterations;
optParam.OutputWeight = outputWeight;

fprintf('\n============================================================\n');
fprintf('STEP 1: Estimate 7 parameters on training run, fixed IC\n');
fprintf('============================================================\n');

modelTrain = setFixedInitialState(model0, x0Train);
modelParam = nlgreyest(zTrain, modelTrain, optParam);
modelParam.Name = 'Final 7-parameter bias model, fixed IC';

fprintf('\nEstimated parameters after training fit:\n');
dispParameters(modelParam);

%% ================================================================
% Optional second refinement
% ================================================================

if runSecondRefinement

    optRefine = optParam;
    optRefine.SearchOptions.MaxIterations = secondRefinementIterations;

    fprintf('\n============================================================\n');
    fprintf('STEP 2: Second parameter refinement, fixed IC\n');
    fprintf('============================================================\n');

    modelParam = nlgreyest(zTrain, modelParam, optRefine);
    modelParam.Name = 'Final 7-parameter bias model, fixed IC, refined';

    fprintf('\nEstimated parameters after second refinement:\n');
    dispParameters(modelParam);
end

%% ================================================================
% Compare training and validation with fixed ICs
% ================================================================

modelTrainFinal = setFixedInitialState(modelParam, x0Train);
modelValFinal   = setFixedInitialState(modelParam, x0Val);

[yHatTrain, fitTrain] = compareUsingModelIC(zTrain, modelTrainFinal);
[yHatVal,   fitVal]   = compareUsingModelIC(zVal,   modelValFinal);

fprintf('\n============================================================\n');
fprintf('FINAL FIT WITH FIXED INITIAL CONDITIONS\n');
fprintf('============================================================\n');

fprintf('\nTraining fit:\n');
dispFitBreakdown(fitTrain, {'theta_1', 'phi_2'});

fprintf('Validation fit:\n');
dispFitBreakdown(fitVal, {'theta_1', 'phi_2'});

fprintf('\nFinal parameter values:\n');
dispParameters(modelParam);

fprintf('\nTraining initial states used:\n');
dispInitialStates(modelTrainFinal);

fprintf('Validation initial states used:\n');
dispInitialStates(modelValFinal);

plotFinalComparison(zTrain, zVal, yHatTrain, yHatVal, tTrain, tVal, fitTrain, fitVal);

fprintf('\nAmplitude diagnostics, training:\n');
printAmplitudeDiagnostics(zTrain, yHatTrain, {'theta_1', 'phi_2'});

fprintf('\nAmplitude diagnostics, validation:\n');
printAmplitudeDiagnostics(zVal, yHatVal, {'theta_1', 'phi_2'});

%% ================================================================
% Optional multistart
% ================================================================

bestModel = modelParam;
bestValMean = meanFitValue(fitVal);

if runMultistart

    fprintf('\n============================================================\n');
    fprintf('OPTIONAL MULTISTART\n');
    fprintf('============================================================\n');

    rng(12);

    optMS = optParam;
    optMS.EstimateCovariance = false;
    optMS.SearchOptions.MaxIterations = multistartIterations;

    for s = 1:nMultistarts

        fprintf('\nMultistart %d / %d\n', s, nMultistarts);

        if s == 1
            modelStart = modelParam;
        else
            modelStart = randomizeFreeParameters(modelParam, 0.35);
        end

        modelStart = setFixedInitialState(modelStart, x0Train);

        try
            modelTry = nlgreyest(zTrain, modelStart, optMS);
            modelTryVal = setFixedInitialState(modelTry, x0Val);

            [~, fitTryVal] = compareUsingModelIC(zVal, modelTryVal);
            meanVal = meanFitValue(fitTryVal);

            fprintf('Validation mean fit: %.3f %%\n', meanVal);
            dispParameters(modelTry);

            if meanVal > bestValMean
                bestValMean = meanVal;
                bestModel = modelTry;
            end

        catch ME
            warning('Multistart %d failed: %s', s, ME.message);
        end
    end
end

%% ================================================================
% Save result
% ================================================================

final7FixedIC.modelParam = modelParam;
final7FixedIC.bestModel = bestModel;
final7FixedIC.modelTrainFinal = modelTrainFinal;
final7FixedIC.modelValFinal = modelValFinal;
final7FixedIC.fitTrain = fitTrain;
final7FixedIC.fitVal = fitVal;
final7FixedIC.trainInfo = trainInfo;
final7FixedIC.valInfo = valInfo;

final7FixedIC.settings.dataFile = dataFile;
final7FixedIC.settings.trainDatasetName = trainDatasetName;
final7FixedIC.settings.validationDatasetName = validationDatasetName;
final7FixedIC.settings.usedStartTime = usedStartTime;
final7FixedIC.settings.usedEndTime = usedEndTime;
final7FixedIC.settings.estimationStartTime = estimationStartTime;
final7FixedIC.settings.inputPolarity = inputPolarity;
final7FixedIC.settings.relativeAngleSign = relativeAngleSign;
final7FixedIC.settings.removeInitialSensorOffsets = removeInitialSensorOffsets;
final7FixedIC.settings.outputWeight = outputWeight;
final7FixedIC.settings.x0Train = x0Train;
final7FixedIC.settings.x0Val = x0Val;

save('final_whitebox_7par_bias_fixed_ic_train_val_result.mat', 'final7FixedIC');

fprintf('\nSaved result to final_whitebox_7par_bias_fixed_ic_train_val_result.mat\n');

%% ================================================================
% Local helper functions
% ================================================================

function [data, info] = loadAndPreprocessRun(raw, datasetName, directNames, ...
    usedStartTime, usedEndTime, inputPolarity, relativeAngleSign, removeInitialSensorOffsets, label)

if strlength(directNames.theta1) > 0
    theta1Name = directNames.theta1;
else
    theta1Name = "theta_1_" + datasetName;
end

if strlength(directNames.theta2) > 0
    theta2Name = directNames.theta2;
else
    theta2Name = "theta_2_" + datasetName;
end

if strlength(directNames.u) > 0
    uName = directNames.u;
else
    uName = "u_" + datasetName;
end

if ~isfield(raw, theta1Name) || ~isfield(raw, theta2Name) || ~isfield(raw, uName)
    fprintf('\nAvailable variables in MAT file:\n');
    disp(string(fieldnames(raw)));
    error('Could not find expected variables for %s run: %s, %s, %s', ...
        label, theta1Name, theta2Name, uName);
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

idxUsed = t_raw >= usedStartTime & t_raw <= usedEndTime;

t      = t_raw(idxUsed);
theta1 = theta1_raw(idxUsed);
phi2   = relativeAngleSign * theta2_raw(idxUsed);
u      = inputPolarity * u_raw(idxUsed);

t = t - t(1);

if removeInitialSensorOffsets
    theta1Offset = theta1(1);
    phi2Offset   = phi2(1);

    theta1 = theta1 - theta1Offset;
    phi2   = phi2   - phi2Offset;
else
    theta1Offset = 0;
    phi2Offset   = 0;
end

theta2_abs = theta1 + phi2;

data.t = t;
data.theta1 = theta1;
data.phi2 = phi2;
data.theta2_abs = theta2_abs;
data.u = u;

info.name = datasetName;
info.theta1Name = theta1Name;
info.theta2Name = theta2Name;
info.uName = uName;
info.theta1Offset = theta1Offset;
info.phi2Offset = phi2Offset;

end

function x0 = fixedInitialStateFromFirstSample(y)

theta1_0 = y(1, 1);
phi2_0   = y(1, 2);
theta2_0 = theta1_0 + phi2_0;

x0 = {theta1_0; theta2_0; 0; 0};

end

function model = configure7ParModel(model)

names = {'p_J1','p_kappa','p_g1','p_g2','p_Ku','p_b1','p_tau0'};
units = {'-','-','1/s^2','1/s^2','torque/input','1/s','torque'};

mins = [0.05; -0.98; 0; 0; 0; 0; -3000];
maxs = [80; 0.98; 800; 800; 20000; 5000; 3000];

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

% Fixed ICs: new run starts with 1 s idle.
model = setinit(model, 'Fixed', {true; true; true; true});

end

function modelOut = setFixedInitialState(modelIn, x0)

modelOut = modelIn;
modelOut = setinit(modelOut, 'Minimum', {-Inf; -Inf; -Inf; -Inf});
modelOut = setinit(modelOut, 'Maximum', { Inf;  Inf;  Inf;  Inf});
modelOut = setinit(modelOut, 'Value', x0);
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

function plotMeasuredRuns(trainData, valData, fs)

figure('Name', 'Measured training and validation data', ...
    'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; plot(trainData.t, trainData.theta1); grid on; ylabel('\theta_1 [rad]'); title('Training \theta_1');
nexttile; plot(valData.t, valData.theta1);   grid on; ylabel('\theta_1 [rad]'); title('Validation \theta_1');
nexttile; plot(trainData.t, trainData.phi2); grid on; ylabel('\phi_2 [rad]'); title('Training \phi_2');
nexttile; plot(valData.t, valData.phi2);     grid on; ylabel('\phi_2 [rad]'); title('Validation \phi_2');
nexttile; plot(trainData.t, trainData.theta2_abs); grid on; ylabel('\theta_{2,abs} [rad]'); title('Training \theta_{2,abs}');
nexttile; plot(valData.t, valData.theta2_abs);     grid on; ylabel('\theta_{2,abs} [rad]'); title('Validation \theta_{2,abs}');
nexttile; plot(trainData.t, trainData.u); grid on; ylabel('u'); xlabel('Time [s]'); title('Training input');
nexttile; plot(valData.t, valData.u);     grid on; ylabel('u'); xlabel('Time [s]'); title('Validation input');

figure('Name', 'Fourier spectra', ...
    'Units', 'normalized', 'Position', [0.08 0.08 0.85 0.75]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
plotSpectrumTile(trainData.theta1, fs, 'Training \theta_1 spectrum');
plotSpectrumTile(valData.theta1,   fs, 'Validation \theta_1 spectrum');
plotSpectrumTile(trainData.phi2,   fs, 'Training \phi_2 spectrum');
plotSpectrumTile(valData.phi2,     fs, 'Validation \phi_2 spectrum');
plotSpectrumTile(trainData.u,      fs, 'Training input spectrum');
plotSpectrumTile(valData.u,        fs, 'Validation input spectrum');

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

figure('Name', 'Final 7-parameter fixed-IC white-box fit', ...
    'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(tTrain, yTrain(:,1), 'DisplayName', 'measured'); hold on;
nTrain1 = min(numel(tTrain), size(yTrainHat, 1));
plot(tTrain(1:nTrain1), yTrainHat(1:nTrain1,1), '--', 'DisplayName', 'model');
grid on; ylabel('\theta_1 [rad]'); title(sprintf('Training \\theta_1, mean %.2f %%', trainFitMean)); legend('Location', 'best');

nexttile;
plot(tTrain, yTrain(:,2), 'DisplayName', 'measured'); hold on;
nTrain2 = min(numel(tTrain), size(yTrainHat, 1));
plot(tTrain(1:nTrain2), yTrainHat(1:nTrain2,2), '--', 'DisplayName', 'model');
grid on; ylabel('\phi_2 [rad]'); title('Training \phi_2'); legend('Location', 'best');

nexttile;
plot(tVal, yVal(:,1), 'DisplayName', 'measured'); hold on;
nVal1 = min(numel(tVal), size(yValHat, 1));
plot(tVal(1:nVal1), yValHat(1:nVal1,1), '--', 'DisplayName', 'model');
grid on; ylabel('\theta_1 [rad]'); xlabel('Time [s]'); title(sprintf('Validation \\theta_1, mean %.2f %%', valFitMean)); legend('Location', 'best');

nexttile;
plot(tVal, yVal(:,2), 'DisplayName', 'measured'); hold on;
nVal2 = min(numel(tVal), size(yValHat, 1));
plot(tVal(1:nVal2), yValHat(1:nVal2,2), '--', 'DisplayName', 'model');
grid on; ylabel('\phi_2 [rad]'); xlabel('Time [s]'); title('Validation \phi_2'); legend('Location', 'best');

sgtitle('Final 7-parameter bias model with fixed initial conditions');

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

function model = randomizeFreeParameters(modelBase, spread)
model = modelBase;
for k = 1:length(model.Parameters)
    if isParameterFixed(model.Parameters(k))
        continue;
    end
    p = double(model.Parameters(k).Value);
    pMin = double(model.Parameters(k).Minimum);
    pMax = double(model.Parameters(k).Maximum);
    if pMin < 0 && pMax > 0
        pNew = p + spread * max(abs(p), 1) * randn;
    elseif p > 0
        pNew = p * exp(spread * randn);
    else
        pNew = pMin + rand * (pMax - pMin);
    end
    pNew = max(pMin, min(pMax, pNew));
    model.Parameters(k).Value = pNew;
end
end
