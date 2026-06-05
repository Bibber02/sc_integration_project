clear;
clear functions;
clc;
close all;

fprintf('\n======================================================\n');
fprintf('Starting passive-link identification: viscous-only workflow, no local helper functions\n');
fprintf('======================================================\n\n');

%% ================================================================
% Configuration - change settings here
% ================================================================

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end

% Change these if your folders are different.
dataFolder   = fullfile(scriptFolder, '..', '..', 'measurement_data');
modelsFolder = fullfile(scriptFolder, '..', 'models');
outputFolder = scriptFolder;

dataFile = fullfile(dataFolder, 'matlab.mat');

% Fallbacks are useful when the script is placed directly in the project root.
if ~exist(dataFile, 'file')
    dataFolder = fullfile(scriptFolder, 'measurement_data');
    dataFile = fullfile(dataFolder, 'matlab.mat');
end
if ~exist(modelsFolder, 'dir')
    modelsFolder = scriptFolder;
end

addpath(modelsFolder);

% idnlgrey expects the callable function name, not the full .m file path.
modelFile = 'passive_rod_viscous_model';

resultMatFile   = fullfile(outputFolder, 'passive_link_viscous_result_v5_no_helpers.mat');
parameterCsvFile = fullfile(outputFolder, 'passive_link_viscous_parameters_v5_no_helpers.csv');
consoleLogFile  = fullfile(outputFolder, 'passive_link_viscous_console_output_v5_no_helpers.txt');

% Plotting. The figure position is intentionally not full-screen.
plotRawData = true;
plotSmoothedPeakData = true;
plotSlicedData = true;
plotFrequencyData = true;
plotFinalFits = true;
figurePosition = [80 80 1050 650];
set(groot, 'defaultFigureWindowStyle', 'normal');
set(groot, 'defaultFigureUnits', 'pixels');
set(groot, 'defaultFigurePosition', figurePosition);

% Smoothing/peak-detection settings for the slicing diagnostic plot only.
smoothingCutoffHz = 5;
smoothingFilterOrder = 3;
peakMinDistanceSeconds = 0.30;

% Measurement definitions.
angles = 15:15:120;
nExp = length(angles);
idxEst = [1 3 5 7];
idxVal = [2 4 6 8];

% Manual slice indices from the original passive-link script.
startIndex = [919 1078 1089 955 1091 1051 851 802];
endIndex   = [1703 2594 3339 3643 3969 4656 4857 4984];

% Fallback starting values. These are overwritten by the data-derived guess
% when the FFT/peak-decay estimate succeeds.
p_b2_start = 0.067;
p_g2_start = 112;

% Initial-state handling. theta_2(0) is fixed; theta_2_dot(0) is estimated.
estimateInitialSpeed = true;
velocityFitSamples = 9;
maxAbsInitialSpeed = 5.0;

% Estimation options.
estimationDisplay = 'Full';
estimationSearchMethod = 'lm';
estimationMaxIterations = 100;
estimateCovariance = true;

%% ================================================================
% Setup
% ================================================================

diary(consoleLogFile);

fprintf('Data file:   %s\n', dataFile);
fprintf('Models path: %s\n', modelsFolder);
fprintf('Model file:  %s\n', modelFile);
fprintf('Result file: %s\n\n', resultMatFile);

whichModel = which(modelFile);
if isempty(whichModel)
    error('MATLAB cannot find %s. Check modelsFolder and addpath(modelsFolder).', modelFile);
else
    fprintf('MATLAB will use model file:\n  %s\n\n', whichModel);
end

%% ================================================================
% Load data
% ================================================================

rawData = load(dataFile);
expData = struct();

for i = 1:nExp
    varName = sprintf('theta_2_%d_degrees', angles(i));
    data = rawData.(varName);

    expData(i).angle_deg = angles(i); %#ok<SAGROW>
    expData(i).t_raw = double(data.Time(:));
    expData(i).theta_raw = double(data.Data(:));
    expData(i).start_index = startIndex(i);
    expData(i).end_index = endIndex(i);
end

%% ================================================================
% Plot raw data
% ================================================================

if plotRawData
    figure('Name', 'Raw passive-link measurements');
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:nExp
        nexttile;
        plot(expData(i).t_raw, expData(i).theta_raw, 'DisplayName', 'Raw data');
        hold on;
        xline(expData(i).t_raw(startIndex(i)), '--', 'Start', 'HandleVisibility', 'off');
        xline(expData(i).t_raw(endIndex(i)), '--', 'End', 'HandleVisibility', 'off');
        title(sprintf('%d degrees', angles(i)));
        ylabel('\theta_2 [rad]');
        xlabel('Time [s]');
    end

    sgtitle('Raw measurements with selected slice limits');
    drawnow;
end

%% ================================================================
% Plot smoothed data and detected peaks for slice checking
% ================================================================

if plotSmoothedPeakData
    figure('Name', 'Smoothed passive-link data with detected peaks');
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:nExp
        TsRaw = expData(i).t_raw(2) - expData(i).t_raw(1);
        fsRaw = 1 / TsRaw;

        [b, a] = butter(smoothingFilterOrder, smoothingCutoffHz/(fsRaw/2), 'low');
        expData(i).theta_smooth = filtfilt(b, a, expData(i).theta_raw);

        minDist = round(peakMinDistanceSeconds * fsRaw);
        [~, locMax] = findpeaks(expData(i).theta_smooth, 'MinPeakDistance', minDist);
        [~, locMin] = findpeaks(-expData(i).theta_smooth, 'MinPeakDistance', minDist);
        peakIndices = sort([locMax; locMin]);
        expData(i).peak_indices = peakIndices;

        nexttile;
        plot(expData(i).t_raw, expData(i).theta_raw, 'DisplayName', 'Raw data');
        hold on;
        plot(expData(i).t_raw, expData(i).theta_smooth, 'DisplayName', 'Smoothed data');
        plot(expData(i).t_raw(peakIndices), expData(i).theta_smooth(peakIndices), 'o', ...
            'DisplayName', 'Detected peaks');
        xline(expData(i).t_raw(startIndex(i)), '--', 'Start', 'HandleVisibility', 'off');
        xline(expData(i).t_raw(endIndex(i)), '--', 'End', 'HandleVisibility', 'off');

        legend('Location', 'best');
        title(sprintf('%d degrees', angles(i)));
        xlabel('Time [s]');
        ylabel('\theta_2 [rad]');
    end

    sgtitle('Smoothed data and detected peaks used to check slice indices');
    drawnow;
end

%% ================================================================
% Slice data and subtract final resting offset
% ================================================================

for i = 1:nExp
    thetaRaw = expData(i).theta_raw;
    tRaw = expData(i).t_raw;

    restSamples = endIndex(i)+1:length(thetaRaw);
    if isempty(restSamples)
        restSamples = max(1, endIndex(i)-100):endIndex(i);
    end
    offsetEstimate = mean(thetaRaw(restSamples));

    expData(i).offsetEstimate = offsetEstimate;
    expData(i).theta = thetaRaw(startIndex(i):endIndex(i)) - offsetEstimate;
    expData(i).t = tRaw(startIndex(i):endIndex(i));
    expData(i).t = expData(i).t - expData(i).t(1);
    expData(i).u = zeros(length(expData(i).theta), 1);

    expData(i).theta0 = expData(i).theta(1);

    if estimateInitialSpeed
        nFit = min(velocityFitSamples, length(expData(i).t));
        tFit = expData(i).t(1:nFit) - expData(i).t(1);
        thetaFit = expData(i).theta(1:nFit);
        if numel(tFit) < 2 || all(abs(tFit) < eps)
            thetaDot0 = 0;
        else
            speedFit = polyfit(tFit, thetaFit, 1);
            thetaDot0 = speedFit(1);
        end
        expData(i).theta_dot0 = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, thetaDot0));
    else
        expData(i).theta_dot0 = 0;
    end
end

Ts = expData(1).t(2) - expData(1).t(1);
fprintf('Sample time Ts = %.6f s\n\n', Ts);

%% ================================================================
% Estimate rough initial parameter guesses from the selected data
% ================================================================

try
    fPeaks = zeros(length(idxEst), 1);
    pbGuesses = zeros(length(idxEst), 1);

    for k = 1:length(idxEst)
        i = idxEst(k);

        y = double(expData(i).theta(:));
        t = double(expData(i).t(:));
        TsLocal = t(2) - t(1);
        fs = 1 / TsLocal;
        N = length(y);

        y0 = y - mean(y);
        Y = fft(y0);
        P2 = abs(Y / N);
        P1 = P2(1:floor(N/2)+1);
        P1(2:end-1) = 2 * P1(2:end-1);
        f = fs * (0:floor(N/2)) / N;

        valid = f > 0.2 & f < 10;
        if ~any(valid)
            error('No valid frequency range found for experiment %d.', i);
        end

        [~, idxMax] = max(P1(valid));
        fValid = f(valid);
        fPeaks(k) = fValid(idxMax);
        omegaD = 2 * pi * fPeaks(k);

        minPeakDist = round(0.35 / TsLocal);
        [pks, ~] = findpeaks(abs(y0), 'MinPeakDistance', minPeakDist);

        if isempty(pks)
            zeta = 0.03;
        else
            pks = pks(pks > 0.02 * max(pks));
            if length(pks) >= 4
                logDec = log(pks(1:end-1) ./ pks(2:end));
                logDec = logDec(isfinite(logDec) & logDec > 0);
                if isempty(logDec)
                    zeta = 0.03;
                else
                    delta = median(logDec);
                    zeta = delta / sqrt(4*pi^2 + delta^2);
                end
            else
                zeta = 0.03;
            end
        end

        pbGuesses(k) = 2 * zeta * omegaD;
    end

    omega0 = 2 * pi * median(fPeaks);
    p_g2_start = max(min(omega0^2, 400), 1);
    p_b2_start = max(min(median(pbGuesses), 30), 0.001);

    fprintf('Data-derived initial guesses:\n');
    fprintf('  p_b2 = %.6g\n', p_b2_start);
    fprintf('  p_g2 = %.6g\n\n', p_g2_start);
catch ME
    warning('Could not compute data-derived initial guesses. Using fallback values instead: %s', ME.message);
    fprintf('Fallback initial guesses:\n');
    fprintf('  p_b2 = %.6g\n', p_b2_start);
    fprintf('  p_g2 = %.6g\n\n', p_g2_start);
end

%% ================================================================
% Plot sliced time-domain data
% ================================================================

if plotSlicedData
    figure('Name', 'Sliced and offset-corrected passive-link data');
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:nExp
        nexttile;
        plot(expData(i).t, expData(i).theta);
        title(sprintf('%d degrees', angles(i)));
        xlabel('Time [s]');
        ylabel('offset-corrected \theta_2 [rad]');
    end

    sgtitle('Sliced and offset-corrected measurements');
    drawnow;
end

%% ================================================================
% Plot sliced frequency-domain data
% ================================================================

if plotFrequencyData
    figure('Name', 'Frequency content of sliced passive-link data');
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:nExp
        y = expData(i).theta(:);
        N = length(y);
        fs = 1 / Ts;

        Y = fft(y - mean(y));
        P2 = abs(Y / N);
        P1 = P2(1:floor(N/2)+1);
        P1(2:end-1) = 2 * P1(2:end-1);
        f = fs * (0:floor(N/2)) / N;

        nexttile;
        plot(f, P1);
        xscale('log');
        title(sprintf('%d degrees', angles(i)));
        xlabel('Frequency [Hz]');
        ylabel('Amplitude');
    end

    sgtitle('Frequency content of sliced measurements');
    drawnow;
end

%% ================================================================
% Build iddata objects
% ================================================================

zEst = [];
zVal = [];
x0Est = [];
x0Val = [];
nEst = 0;
nVal = 0;

for i = 1:nExp
    z = iddata(expData(i).theta(:), expData(i).u(:), Ts);
    z.Name = sprintf('passive_%03d_deg', angles(i));
    z.InputName = {'u'};
    z.InputUnit = {'-'};
    z.OutputName = {'theta_2'};
    z.OutputUnit = {'rad'};
    z.TimeUnit = 's';

    x0 = [expData(i).theta0; expData(i).theta_dot0];

    if ismember(i, idxEst)
        nEst = nEst + 1;
        if isempty(zEst)
            zEst = z;
        else
            zEst = merge(zEst, z);
        end
        x0Est(:, nEst) = x0; %#ok<SAGROW>
    else
        nVal = nVal + 1;
        if isempty(zVal)
            zVal = z;
        else
            zVal = merge(zVal, z);
        end
        x0Val(:, nVal) = x0; %#ok<SAGROW>
    end
end

fprintf('Number of estimation experiments: %d\n', nEst);
fprintf('Number of validation experiments: %d\n\n', nVal);

%% ================================================================
% Grey-box model setup
% ================================================================

order = [1 1 2];
parameters = {p_b2_start; p_g2_start};
initialStates = {x0Est(1, :); x0Est(2, :)};
model0 = idnlgrey(modelFile, order, parameters, initialStates, 0);

model0.TimeUnit = 's';
model0.InputName = {'u'};
model0.InputUnit = {'-'};
model0.OutputName = {'theta_2'};
model0.OutputUnit = {'rad'};

model0.Parameters(1).Name = 'p_b2';
model0.Parameters(1).Unit = '1/s';
model0.Parameters(1).Minimum = 0;
model0.Parameters(1).Maximum = 30;
model0.Parameters(1).Fixed = false;

model0.Parameters(2).Name = 'p_g2';
model0.Parameters(2).Unit = '1/s^2';
model0.Parameters(2).Minimum = 1;
model0.Parameters(2).Maximum = 400;
model0.Parameters(2).Fixed = false;

model0.InitialStates(1).Name = 'theta_2';
model0.InitialStates(1).Unit = 'rad';
model0.InitialStates(2).Name = 'theta_2_dot';
model0.InitialStates(2).Unit = 'rad/s';

model0 = setinit(model0, 'Fixed', {true(1, nEst); false(1, nEst)});

%% ================================================================
% Estimate model
% ================================================================

compareOpt = compareOptions;
compareOpt.InitialCondition = 'estimate';

opt = nlgreyestOptions;
opt.Display = estimationDisplay;
opt.EstimateCovariance = estimateCovariance;
opt.SearchMethod = estimationSearchMethod;
opt.SearchOptions.MaxIterations = estimationMaxIterations;

fprintf('\n======================================================\n');
fprintf('Estimating viscous passive-link model\n');
fprintf('======================================================\n\n');

model_viscous = nlgreyest(zEst, model0, opt);
model_viscous.Name = 'Final passive-link viscous model';

fprintf('\nFinal viscous parameter values:\n');
for k = 1:length(model_viscous.Parameters)
    fprintf('%-12s = %12.6g\n', model_viscous.Parameters(k).Name, model_viscous.Parameters(k).Value);
end

%% ================================================================
% Validation model with validation initial states
% ================================================================

parametersVal = cell(length(model_viscous.Parameters), 1);
for k = 1:length(model_viscous.Parameters)
    parametersVal{k} = model_viscous.Parameters(k).Value;
end

model_viscous_val = idnlgrey(modelFile, order, parametersVal, {x0Val(1, :); x0Val(2, :)}, 0);
model_viscous_val.Name = model_viscous.Name;
model_viscous_val.TimeUnit = model_viscous.TimeUnit;
model_viscous_val.InputName = model_viscous.InputName;
model_viscous_val.InputUnit = model_viscous.InputUnit;
model_viscous_val.OutputName = model_viscous.OutputName;
model_viscous_val.OutputUnit = model_viscous.OutputUnit;

for k = 1:length(model_viscous.Parameters)
    model_viscous_val.Parameters(k).Name = model_viscous.Parameters(k).Name;
    model_viscous_val.Parameters(k).Unit = model_viscous.Parameters(k).Unit;
    model_viscous_val.Parameters(k).Minimum = model_viscous.Parameters(k).Minimum;
    model_viscous_val.Parameters(k).Maximum = model_viscous.Parameters(k).Maximum;
    model_viscous_val.Parameters(k).Fixed = model_viscous.Parameters(k).Fixed;
end

model_viscous_val.InitialStates(1).Name = model_viscous.InitialStates(1).Name;
model_viscous_val.InitialStates(1).Unit = model_viscous.InitialStates(1).Unit;
model_viscous_val.InitialStates(2).Name = model_viscous.InitialStates(2).Name;
model_viscous_val.InitialStates(2).Unit = model_viscous.InitialStates(2).Unit;
model_viscous_val = setinit(model_viscous_val, 'Fixed', {true(1, nVal); false(1, nVal)});

%% ================================================================
% Validation and fit statistics
% ================================================================

[yEst, fitEst] = compare(zEst, model_viscous, compareOpt);
[yVal, fitVal] = compare(zVal, model_viscous_val, compareOpt);

fitEstValues = [];
if isnumeric(fitEst)
    fitEstValues = fitEst(:);
elseif iscell(fitEst)
    for k = 1:numel(fitEst)
        if isnumeric(fitEst{k})
            fitEstValues = [fitEstValues; fitEst{k}(:)]; %#ok<AGROW>
        elseif iscell(fitEst{k})
            for kk = 1:numel(fitEst{k})
                if isnumeric(fitEst{k}{kk})
                    fitEstValues = [fitEstValues; fitEst{k}{kk}(:)]; %#ok<AGROW>
                end
            end
        end
    end
end
fitEstValues = fitEstValues(isfinite(fitEstValues));
if isempty(fitEstValues)
    meanFitEst = NaN;
else
    meanFitEst = mean(fitEstValues);
end

fitValValues = [];
if isnumeric(fitVal)
    fitValValues = fitVal(:);
elseif iscell(fitVal)
    for k = 1:numel(fitVal)
        if isnumeric(fitVal{k})
            fitValValues = [fitValValues; fitVal{k}(:)]; %#ok<AGROW>
        elseif iscell(fitVal{k})
            for kk = 1:numel(fitVal{k})
                if isnumeric(fitVal{k}{kk})
                    fitValValues = [fitValValues; fitVal{k}{kk}(:)]; %#ok<AGROW>
                end
            end
        end
    end
end
fitValValues = fitValValues(isfinite(fitValValues));
if isempty(fitValValues)
    meanFitVal = NaN;
else
    meanFitVal = mean(fitValValues);
end

fprintf('\nViscous fit, per estimation experiment:\n');
disp(fitEst);
fprintf('\nViscous fit, per validation experiment:\n');
disp(fitVal);
fprintf('\nMean estimation fit: %.2f %%\n', meanFitEst);
fprintf('Mean validation fit: %.2f %%\n', meanFitVal);

%% ================================================================
% Parameter uncertainty table
% ================================================================

nPar = length(model_viscous.Parameters);
parameter = cell(nPar, 1);
value = zeros(nPar, 1);
fixed = false(nPar, 1);
variance = NaN(nPar, 1);
stdDev = NaN(nPar, 1);
relativeStdPercent = NaN(nPar, 1);

for k = 1:nPar
    parameter{k} = model_viscous.Parameters(k).Name;
    value(k) = model_viscous.Parameters(k).Value;
    fixed(k) = model_viscous.Parameters(k).Fixed;
end

freeIdx = find(~fixed);
try
    covFree = getcov(model_viscous, 'value', 'free');
    for k = 1:length(freeIdx)
        variance(freeIdx(k)) = covFree(k, k);
        stdDev(freeIdx(k)) = sqrt(max(covFree(k, k), 0));
        if abs(value(freeIdx(k))) > eps
            relativeStdPercent(freeIdx(k)) = 100 * stdDev(freeIdx(k)) / abs(value(freeIdx(k)));
        end
    end
catch ME
    warning('Could not compute covariance table: %s', ME.message);
end

resultTable = table(parameter, value, fixed, variance, stdDev, relativeStdPercent, ...
    'VariableNames', {'Parameter', 'Value', 'Fixed', 'Variance', 'StdDev', 'RelativeStdPercent'});

fprintf('\nFinal parameter table:\n');
disp(resultTable);

%% ================================================================
% Convert compare outputs to cell arrays for plotting
% ================================================================

yEstCell = cell(numel(idxEst), 1);
if iscell(yEst)
    for k = 1:numel(yEst)
        yRaw = get(yEst{k}, 'OutputData');
        if iscell(yRaw)
            yRaw = yRaw{1};
        end
        yEstCell{k} = double(yRaw(:));
    end
else
    yRaw = get(yEst, 'OutputData');
    if iscell(yRaw)
        for k = 1:numel(yRaw)
            yEstCell{k} = double(yRaw{k}(:));
        end
    else
        yEstCell{1} = double(yRaw(:));
    end
end

yValCell = cell(numel(idxVal), 1);
if iscell(yVal)
    for k = 1:numel(yVal)
        yRaw = get(yVal{k}, 'OutputData');
        if iscell(yRaw)
            yRaw = yRaw{1};
        end
        yValCell{k} = double(yRaw(:));
    end
else
    yRaw = get(yVal, 'OutputData');
    if iscell(yRaw)
        for k = 1:numel(yRaw)
            yValCell{k} = double(yRaw{k}(:));
        end
    else
        yValCell{1} = double(yRaw(:));
    end
end

%% ================================================================
% Final plots
% ================================================================

if plotFinalFits
    figure('Name', 'Viscous identification fit');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:length(idxEst)
        i = idxEst(k);
        nexttile;
        measured = expData(i).theta(:);
        if k <= numel(yEstCell) && ~isempty(yEstCell{k})
            modelled = yEstCell{k}(:);
        else
            modelled = NaN(size(measured));
        end
        n = min(length(measured), length(modelled));
        plot(expData(i).t(1:n), measured(1:n), 'DisplayName', 'measured');
        hold on;
        plot(expData(i).t(1:n), modelled(1:n), '--', 'DisplayName', 'model');
        if k <= length(fitEstValues) && isfinite(fitEstValues(k))
            fitText = sprintf(', fit %.1f%%', fitEstValues(k));
        else
            fitText = '';
        end
        title(sprintf('%d deg%s', angles(i), fitText));
        xlabel('Time [s]');
        ylabel('\theta_2 [rad]');
        legend('Location', 'best');
    end
    sgtitle('Viscous model: estimation data');

    figure('Name', 'Viscous validation fit');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:length(idxVal)
        i = idxVal(k);
        nexttile;
        measured = expData(i).theta(:);
        if k <= numel(yValCell) && ~isempty(yValCell{k})
            modelled = yValCell{k}(:);
        else
            modelled = NaN(size(measured));
        end
        n = min(length(measured), length(modelled));
        plot(expData(i).t(1:n), measured(1:n), 'DisplayName', 'measured');
        hold on;
        plot(expData(i).t(1:n), modelled(1:n), '--', 'DisplayName', 'model');
        if k <= length(fitValValues) && isfinite(fitValValues(k))
            fitText = sprintf(', fit %.1f%%', fitValValues(k));
        else
            fitText = '';
        end
        title(sprintf('%d deg%s', angles(i), fitText));
        xlabel('Time [s]');
        ylabel('\theta_2 [rad]');
        legend('Location', 'best');
    end
    sgtitle('Viscous model: validation data');
    drawnow;
end

%% ================================================================
% Save results
% ================================================================

passiveID.model_viscous = model_viscous;
passiveID.model_viscous_val = model_viscous_val;
passiveID.fitEst = fitEst;
passiveID.fitVal = fitVal;
passiveID.meanFitEst = meanFitEst;
passiveID.meanFitVal = meanFitVal;
passiveID.resultTable = resultTable;
passiveID.expData = expData;
passiveID.idxEst = idxEst;
passiveID.idxVal = idxVal;
passiveID.angles = angles;
passiveID.config.dataFile = dataFile;
passiveID.config.modelsFolder = modelsFolder;
passiveID.config.startIndex = startIndex;
passiveID.config.endIndex = endIndex;
passiveID.config.figurePosition = figurePosition;

save(resultMatFile, 'passiveID');
writetable(resultTable, parameterCsvFile);

fprintf('\nSaved result file: %s\n', resultMatFile);
fprintf('Saved parameter table: %s\n', parameterCsvFile);

diary off;
