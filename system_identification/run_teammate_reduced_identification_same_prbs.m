clear;
clear functions;
clc;
close all;

fprintf('Running run_teammate_reduced_identification_same_prbs.m\n');

%% ================================================================
% User settings: SAME DATASET/SPLIT AS run_compact_relative_identification.m
% ================================================================

dataFile = "identification_data\whole_system_identification\prbs_long.mat";
datasetName = "PRBS_Long2_Amp0_15";

usedEndTime = 90.0;
trainStart  = 0.0;
trainEnd    = 60.0;
valStart    = 60.0;
valEnd      = 90.0;

% Replace this with the measured first-rod length used in your report.
known_l1 = 0.300;       % [m]
known_g  = 9.81;        % [m/s^2]

inputPolarity = -1;
relativeAngleSign = 1;
removeInitialSensorOffsets = true;

% false = 6 estimated parameters: p_J, p_g1, p_g2, p_b1, p_b2, p_u
% true  = 7 estimated parameters: also estimate p_0 = tau_0/beta
estimateTorqueOffset = false;

outputWeight = diag([1, 8]);

parameterMaxIterations = 120;
icOnlyMaxIterations    = 50;

% More conservative than the earlier script. The measured positions are
% fixed to the first sample of each segment; only velocities are estimated.
fixInitialPositions = true;
velocityLimit = 15.0;       % rad/s

%% ================================================================
% Load and preprocess PRBS data
% ================================================================

raw = load(dataFile);

theta1Name = char("theta_1_" + datasetName);
theta2Name = char("theta_2_" + datasetName);  % relative passive sensor
uName      = char("u_"       + datasetName);

requiredNames = {theta1Name, theta2Name, uName};
for k = 1:numel(requiredNames)
    if ~isfield(raw, requiredNames{k})
        error('Missing variable in data file: %s', requiredNames{k});
    end
end

[t_raw, theta1_raw] = timeseriesToVectors(raw.(theta1Name));
[~,     theta2_raw] = timeseriesToVectors(raw.(theta2Name));
[~,     u_raw]      = timeseriesToVectors(raw.(uName));

theta2_raw = relativeAngleSign * theta2_raw;
u_raw      = inputPolarity * u_raw;

N = min([numel(t_raw), numel(theta1_raw), numel(theta2_raw), numel(u_raw)]);
t_raw      = t_raw(1:N);
theta1_raw = theta1_raw(1:N);
theta2_raw = theta2_raw(1:N);
u_raw      = u_raw(1:N);

t_raw = t_raw - t_raw(1);
idxUsed = t_raw <= usedEndTime;

t      = t_raw(idxUsed);
theta1 = theta1_raw(idxUsed);
theta2 = theta2_raw(idxUsed);
u      = u_raw(idxUsed);

if removeInitialSensorOffsets
    theta1Offset = theta1(1);
    theta2Offset = theta2(1);
    theta1 = theta1 - theta1Offset;
    theta2 = theta2 - theta2Offset;
else
    theta1Offset = 0;
    theta2Offset = 0;
end

t = t - t(1);
Ts = mean(diff(t));

theta2AbsForPlot = theta1 + theta2;

fprintf('\nSelected dataset: %s\n', datasetName);
fprintf('Sampling time Ts: %.8f s\n', Ts);
fprintf('Used duration: %.3f s\n', t(end));
fprintf('Input polarity: %+d\n', inputPolarity);
fprintf('Relative angle sign: %+d\n', relativeAngleSign);
fprintf('Known l1: %.8g m\n', known_l1);
fprintf('Known g: %.8g m/s^2\n', known_g);
fprintf('Estimate torque offset p_0: %d\n', estimateTorqueOffset);
fprintf('Removed theta_1 offset: %.8g rad\n', theta1Offset);
fprintf('Removed theta_2 offset: %.8g rad\n', theta2Offset);
fprintf('Offsets are passed into the grey-box model for sin/cos terms.\n');

plotMeasuredData(t, theta1, theta2, theta2AbsForPlot, u);

%% ================================================================
% Build training and validation iddata objects
% ================================================================

idxTrain = t >= trainStart & t < trainEnd;
idxVal   = t >= valStart   & t <= valEnd;

tTrain = t(idxTrain) - t(find(idxTrain, 1, 'first'));
tVal   = t(idxVal)   - t(find(idxVal,   1, 'first'));

yTrain = [theta1(idxTrain), theta2(idxTrain)];
uTrain = u(idxTrain);

yVal = [theta1(idxVal), theta2(idxVal)];
uVal = u(idxVal);

zTrain = makeIdData(yTrain, uTrain, Ts);
zVal   = makeIdData(yVal,   uVal,   Ts);

fprintf('\nTraining samples:   %d, duration %.3f s\n', numel(uTrain), tTrain(end));
fprintf('Validation samples: %d, duration %.3f s\n', numel(uVal), tVal(end));

%% ================================================================
% Create teammate reduced grey-box model
% ================================================================

modelFile = 'teammate_reduced_pendulum_model_offset';
Order = [2 1 4];
TsModel = 0;

% Initial guesses. p_J is p_alpha - p_c^2.
p0.p_g2 = 112.0;
p0.p_c  = (known_l1 / known_g) * p0.p_g2;
p0.p_J  = 8.0;
p0.p_g1 = 100.0;
p0.p_b1 = 650.0;
p0.p_b2 = 0.10;
p0.p_u  = 4000.0;
p0.p_0  = 0.0;

omega1_0 = estimateInitialVelocity(yTrain(:,1), Ts, velocityLimit);
omega2_0 = estimateInitialVelocity(yTrain(:,2), Ts, velocityLimit);

InitialStates = {yTrain(1,1); yTrain(1,2); omega1_0; omega2_0};
Parameters = {
    p0.p_J
    p0.p_g1
    p0.p_g2
    p0.p_b1
    p0.p_b2
    p0.p_u
    p0.p_0
    known_l1
    known_g
    theta1Offset
    theta2Offset
};

model0 = idnlgrey(modelFile, Order, Parameters, InitialStates, TsModel);
model0 = configureTeammateReducedModel(model0, estimateTorqueOffset);
model0 = configureInitialStateBounds(model0, yTrain, Ts, velocityLimit, fixInitialPositions, 'training');

fprintf('\nInitial model parameters:\n');
dispParameters(model0);
printRecoveredParameters(model0);
fprintf('Initial states:\n');
dispInitialStates(model0);

%% ================================================================
% Step 1: estimate parameters and training initial velocities
% ================================================================

optParam = nlgreyestOptions;
optParam.Display = 'on';
optParam.EstimateCovariance = true;
optParam.SearchOptions.MaxIterations = parameterMaxIterations;
optParam.OutputWeight = outputWeight;

fprintf('\n============================================================\n');
fprintf('STEP 1: Estimate teammate reduced parameters and training ICs\n');
fprintf('============================================================\n');

modelParam = nlgreyest(zTrain, model0, optParam);
modelParam.Name = 'Teammate reduced relative-coordinate model, same PRBS dataset';

fprintf('\nEstimated parameters:\n');
dispParameters(modelParam);
printRecoveredParameters(modelParam);
fprintf('Estimated training initial states:\n');
dispInitialStates(modelParam);

%% ================================================================
% Step 2: re-estimate ICs only for train and validation
% ================================================================

optIC = nlgreyestOptions;
optIC.Display = 'on';
optIC.EstimateCovariance = false;
optIC.SearchOptions.MaxIterations = icOnlyMaxIterations;
optIC.OutputWeight = outputWeight;

fprintf('\n============================================================\n');
fprintf('STEP 2: Re-estimate training ICs only, parameters fixed\n');
fprintf('============================================================\n');
[modelTrainIC, fitTrain] = estimateICOnlyAndCompare( ...
    zTrain, modelParam, yTrain, Ts, velocityLimit, fixInitialPositions, 'training', optIC);

fprintf('\n============================================================\n');
fprintf('STEP 3: Estimate validation ICs only, parameters fixed\n');
fprintf('============================================================\n');
[modelValIC, fitVal] = estimateICOnlyAndCompare( ...
    zVal, modelParam, yVal, Ts, velocityLimit, fixInitialPositions, 'validation', optIC);

[yHatTrain, fitTrain] = compareUsingModelIC(zTrain, modelTrainIC);
[yHatVal,   fitVal]   = compareUsingModelIC(zVal,   modelValIC);

fprintf('\nFinal training fit:\n');
dispFitBreakdown(fitTrain, {'theta_1', 'theta_2'});
fprintf('\nFinal validation fit:\n');
dispFitBreakdown(fitVal, {'theta_1', 'theta_2'});

fprintf('\nFinal parameter values:\n');
dispParameters(modelParam);
printRecoveredParameters(modelParam);

plotFitComparison(zTrain, zVal, yHatTrain, yHatVal, tTrain, tVal, fitTrain, fitVal, ...
    'Teammate reduced model on same PRBS dataset');

%% ================================================================
% Save result
% ================================================================

teammateReducedSamePRBS.modelParam = modelParam;
teammateReducedSamePRBS.modelTrainIC = modelTrainIC;
teammateReducedSamePRBS.modelValIC = modelValIC;
teammateReducedSamePRBS.fitTrain = fitTrain;
teammateReducedSamePRBS.fitVal = fitVal;
teammateReducedSamePRBS.settings.dataFile = dataFile;
teammateReducedSamePRBS.settings.datasetName = datasetName;
teammateReducedSamePRBS.settings.known_l1 = known_l1;
teammateReducedSamePRBS.settings.known_g = known_g;
teammateReducedSamePRBS.settings.theta1Offset = theta1Offset;
teammateReducedSamePRBS.settings.theta2Offset = theta2Offset;
teammateReducedSamePRBS.settings.estimateTorqueOffset = estimateTorqueOffset;
teammateReducedSamePRBS.settings.fixInitialPositions = fixInitialPositions;
teammateReducedSamePRBS.settings.velocityLimit = velocityLimit;
teammateReducedSamePRBS.settings.outputWeight = outputWeight;

save('teammate_reduced_same_prbs_result.mat', 'teammateReducedSamePRBS');
fprintf('\nSaved result to teammate_reduced_same_prbs_result.mat\n');

%% ================================================================
% Local helper functions
% ================================================================

function z = makeIdData(y, u, Ts)
z = iddata(y, u, Ts);
z.InputName  = {'u'};
z.OutputName = {'theta_1', 'theta_2'};
z.InputUnit  = {'input'};
z.OutputUnit = {'rad', 'rad'};
z.TimeUnit   = 's';
end

function model = configureTeammateReducedModel(model, estimateTorqueOffset)
names = {'p_J', 'p_g1', 'p_g2', 'p_b1', 'p_b2', 'p_u', 'p_0', 'l1', 'g', 'theta1_offset', 'theta2_offset'};
units = {'-', '1/s^2', '1/s^2', '1/s', '1/s', 'torque/input/beta', 'torque/beta', 'm', 'm/s^2', 'rad', 'rad'};
mins  = [0.001, 0, 0, 0, 0, 0, -5000, model.Parameters(8).Value, model.Parameters(9).Value, model.Parameters(10).Value, model.Parameters(11).Value];
maxs  = [200.0, 800, 800, 6000, 50, 20000, 5000, model.Parameters(8).Value, model.Parameters(9).Value, model.Parameters(10).Value, model.Parameters(11).Value];
fixed = [false, false, false, false, false, false, ~estimateTorqueOffset, true, true, true, true];
for k = 1:numel(names)
    model.Parameters(k).Name = names{k};
    model.Parameters(k).Unit = units{k};
    model.Parameters(k).Minimum = mins(k);
    model.Parameters(k).Maximum = maxs(k);
    model.Parameters(k).Fixed = fixed(k);
end
model.TimeUnit = 's';
model.InputName = {'u'};
model.InputUnit = {'input'};
model.OutputName = {'theta_1', 'theta_2'};
model.OutputUnit = {'rad', 'rad'};
model.InitialStates(1).Name = 'theta_1';
model.InitialStates(1).Unit = 'rad';
model.InitialStates(2).Name = 'theta_2';
model.InitialStates(2).Unit = 'rad';
model.InitialStates(3).Name = 'omega_1';
model.InitialStates(3).Unit = 'rad/s';
model.InitialStates(4).Name = 'omega_2';
model.InitialStates(4).Unit = 'rad/s';
end

function model = configureInitialStateBounds(model, y, Ts, velocityLimit, fixPositions, label)
theta1_0 = y(1, 1);
theta2_0 = y(1, 2);
omega1_0 = estimateInitialVelocity(y(:,1), Ts, velocityLimit);
omega2_0 = estimateInitialVelocity(y(:,2), Ts, velocityLimit);
values = {theta1_0; theta2_0; omega1_0; omega2_0};
if fixPositions
    posTol = 1e-9;
    mins  = {theta1_0-posTol; theta2_0-posTol; -velocityLimit; -velocityLimit};
    maxs  = {theta1_0+posTol; theta2_0+posTol;  velocityLimit;  velocityLimit};
    fixed = {true; true; false; false};
else
    positionMargin = 0.25;
    mins  = {theta1_0-positionMargin; theta2_0-positionMargin; -velocityLimit; -velocityLimit};
    maxs  = {theta1_0+positionMargin; theta2_0+positionMargin;  velocityLimit;  velocityLimit};
    fixed = {false; false; false; false};
end
model = setinit(model, 'Minimum', {-Inf; -Inf; -Inf; -Inf});
model = setinit(model, 'Maximum', { Inf;  Inf;  Inf;  Inf});
model = setinit(model, 'Value', values);
model = setinit(model, 'Minimum', mins);
model = setinit(model, 'Maximum', maxs);
model = setinit(model, 'Fixed', fixed);
fprintf('Configured %s IC bounds. Positions fixed: %d.\n', label, fixPositions);
end

function [modelIC, fit] = estimateICOnlyAndCompare(z, modelParam, y, Ts, velocityLimit, fixPositions, label, optIC)
modelIC = modelParam;
for k = 1:length(modelIC.Parameters)
    modelIC.Parameters(k).Fixed = true;
end
modelIC = configureInitialStateBounds(modelIC, y, Ts, velocityLimit, fixPositions, label);
modelIC = nlgreyest(z, modelIC, optIC);
modelIC.Name = sprintf('%s IC-only teammate reduced model', label);
[~, fit] = compareUsingModelIC(z, modelIC);
end

function v0 = estimateInitialVelocity(y, Ts, limit)
y = y(:);
n = min(10, numel(y));
if n < 2 || Ts <= 0
    v0 = 0;
else
    tt = (0:n-1)' * Ts;
    p = polyfit(tt, y(1:n), 1);
    v0 = p(1);
end
v0 = max(min(v0, limit), -limit);
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

function [t, y] = timeseriesToVectors(signal)
if isa(signal, 'timeseries')
    t = double(signal.Time(:));
    y = double(signal.Data(:));
elseif isstruct(signal) && isfield(signal, 'Time') && isfield(signal, 'Data')
    t = double(signal.Time(:));
    y = double(signal.Data(:));
else
    error('Expected a timeseries object or a struct with Time and Data fields.');
end
end

function plotMeasuredData(t, theta1, theta2, theta2Abs, u)
figure('Name', 'Measured PRBS data', 'Units', 'normalized', 'Position', [0.06 0.06 0.88 0.82]);
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; plot(t, theta1); grid on; ylabel('\theta_1 [rad]'); title('Measured rod-1 local angle');
nexttile; plot(t, theta2); grid on; ylabel('\theta_2 [rad]'); title('Measured relative passive angle');
nexttile; plot(t, theta2Abs); grid on; ylabel('\theta_1+\theta_2 [rad]'); title('Local absolute rod-2 angle for inspection');
nexttile; plot(t, u); grid on; ylabel('u'); xlabel('Time [s]'); title('Input');
end

function plotFitComparison(zTrain, zVal, yHatTrain, yHatVal, tTrain, tVal, fitTrain, fitVal, figTitle)
yTrain = getOutputDataMatrix(zTrain); yVal = getOutputDataMatrix(zVal);
yTrainHat = getOutputDataMatrix(yHatTrain); yValHat = getOutputDataMatrix(yHatVal);
trainMean = meanFitValue(fitTrain); valMean = meanFitValue(fitVal);
figure('Name', figTitle, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; plot(tTrain, yTrain(:,1), 'DisplayName', 'measured'); hold on; n=min(numel(tTrain),size(yTrainHat,1)); plot(tTrain(1:n), yTrainHat(1:n,1), '--', 'DisplayName', 'model'); grid on; ylabel('\theta_1 [rad]'); title(sprintf('Training \theta_1, mean %.2f %%', trainMean)); legend('Location','best');
nexttile; plot(tTrain, yTrain(:,2), 'DisplayName', 'measured'); hold on; n=min(numel(tTrain),size(yTrainHat,1)); plot(tTrain(1:n), yTrainHat(1:n,2), '--', 'DisplayName', 'model'); grid on; ylabel('\theta_2 [rad]'); title('Training \theta_2'); legend('Location','best');
nexttile; plot(tVal, yVal(:,1), 'DisplayName', 'measured'); hold on; n=min(numel(tVal),size(yValHat,1)); plot(tVal(1:n), yValHat(1:n,1), '--', 'DisplayName', 'model'); grid on; ylabel('\theta_1 [rad]'); xlabel('Time [s]'); title(sprintf('Validation \theta_1, mean %.2f %%', valMean)); legend('Location','best');
nexttile; plot(tVal, yVal(:,2), 'DisplayName', 'measured'); hold on; n=min(numel(tVal),size(yValHat,1)); plot(tVal(1:n), yValHat(1:n,2), '--', 'DisplayName', 'model'); grid on; ylabel('\theta_2 [rad]'); xlabel('Time [s]'); title('Validation \theta_2'); legend('Location','best');
sgtitle(figTitle);
end

function y = getOutputDataMatrix(dataObject)
try
    y = get(dataObject, 'OutputData');
catch
    y = dataObject.OutputData;
end
if iscell(y), y = y{1}; end
y = double(y); if isvector(y), y = y(:); end
end

function dispParameters(model)
for k = 1:length(model.Parameters)
    fprintf('  %-14s = %.8g', model.Parameters(k).Name, model.Parameters(k).Value);
    if isFixed(model.Parameters(k).Fixed), fprintf('  fixed\n'); else, fprintf('  free\n'); end
end
end

function printRecoveredParameters(model)
p_J = getPar(model, 'p_J');
p_g2 = getPar(model, 'p_g2');
l1 = getPar(model, 'l1');
g = getPar(model, 'g');
p_c = (l1/g) * p_g2;
p_alpha = p_J + p_c^2;
fprintf('  recovered p_c     = %.8g  fixed by p_g2*l1/g\n', p_c);
fprintf('  recovered p_alpha = %.8g  equal to p_J + p_c^2\n', p_alpha);
end

function value = getPar(model, name)
for k = 1:length(model.Parameters)
    if strcmp(model.Parameters(k).Name, name)
        value = double(model.Parameters(k).Value);
        return;
    end
end
error('Parameter %s not found.', name);
end

function dispInitialStates(model)
for k = 1:length(model.InitialStates)
    fprintf('  %-8s = %.8g', model.InitialStates(k).Name, model.InitialStates(k).Value);
    if isFixed(model.InitialStates(k).Fixed), fprintf('  fixed\n'); else, fprintf('  free\n'); end
end
end

function tf = isFixed(value)
if iscell(value), value = value{1}; end
tf = all(value(:));
end

function dispFitBreakdown(fit, outputNames)
values = collectFitValues(fit); values = values(isfinite(values));
for k = 1:min(numel(values), numel(outputNames))
    fprintf('  %-8s = %.3f %%\n', outputNames{k}, values(k));
end
if isempty(values), fprintf('  mean     = NaN %%\n'); else, fprintf('  mean     = %.3f %%\n', mean(values)); end
end

function m = meanFitValue(fit)
values = collectFitValues(fit); values = values(isfinite(values));
if isempty(values), m = NaN; else, m = mean(values); end
end

function values = collectFitValues(fit)
if isnumeric(fit)
    values = fit(:);
elseif iscell(fit)
    values = [];
    for k = 1:numel(fit), values = [values; collectFitValues(fit{k})]; %#ok<AGROW>
    end
else
    values = NaN;
end
end
