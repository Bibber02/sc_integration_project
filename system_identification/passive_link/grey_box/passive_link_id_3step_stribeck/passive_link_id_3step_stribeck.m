clear;
clear functions;
clc;
close all;

fprintf('\n======================================================\n');
fprintf('Starting passive-link identification: 3-step Stribeck workflow, no local helper functions\n');
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

% idnlgrey expects callable function names, not full .m file paths.
viscousModelFile  = 'passive_rod_viscous_model';
coulombModelFile  = 'passive_rod_coulomb_model';
stribeckModelFile = 'passive_rod_stribeck_report_model';

resultMatFile = fullfile(outputFolder, 'passive_link_3step_stribeck_result_v5_no_helpers.mat');
uncertaintyCsvFile = fullfile(outputFolder, 'passive_link_3step_stribeck_parameter_uncertainty_v5_no_helpers.csv');
identifiabilityCsvFile = fullfile(outputFolder, 'passive_link_3step_stribeck_identifiability_v5_no_helpers.csv');
consoleLogFile = fullfile(outputFolder, 'passive_link_3step_stribeck_console_output_v5_no_helpers.txt');

% Plotting. The figure position is intentionally not full-screen.
plotRawData = true;
plotSmoothedPeakData = true;
plotSlicedData = true;
plotFrequencyData = true;
plotTailSelection = true;
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
idxAll = 1:nExp;

% Manual slice indices from the original passive-link script.
startIndex = [919 1078 1089 955 1091 1051 851 802];
endIndex   = [1703 2594 3339 3643 3969 4656 4857 4984];

% Tail selection for low-speed Stribeck fitting.
thetaTailMax = 0.30;             % rad
minTailDuration = 3.0;           % s
fallbackTailFraction = 0.40;     % last 40 percent if amplitude logic gives too little data
envelopeWindowSeconds = 0.60;    % s moving maximum window

% Runtime and model-selection settings.
runCoulombMultistart = false;
nStartsCoulomb = 1;
stribeckVsCandidates = [2.50 2.75 3.00 3.25 3.50];
nStartsStribeckPerVs = 1;
tailValidationWeight = 0.30;
fullValidationWeight = 0.70;

% Initial-state handling. theta_2(0) is fixed; theta_2_dot(0) is estimated.
estimateInitialSpeed = true;
velocityFitSamples = 9;
maxAbsInitialSpeed = 5.0;

% Estimation options.
estimationDisplay = 'on';
multistartDisplay = 'off';
maxIterationsSimple = 60;
maxIterationsCoulombMultistart = 50;
maxIterationsStribeckGrid = 90;
estimateCovariance = true;

% Parameter metadata.
viscousParNames = {'p_b2', 'p_g2'};
viscousParUnits = {'1/s', '1/s^2'};
viscousParMin = [0, 1];
viscousParMax = [30, 400];

coulombParNames = {'p_b2', 'p_g2', 'p_c2', 'eps_v2'};
coulombParUnits = {'1/s', '1/s^2', 'rad/s^2', 'rad/s'};
coulombParMin = [0, 1, 0, 0.002];
coulombParMax = [30, 400, 80, 0.20];

stribeckParNames = {'p_b2', 'p_g2', 'p_c2', 'p_sdelta2', 'v_s2', 'eps_v2'};
stribeckParUnits = {'1/s', '1/s^2', 'rad/s^2', 'rad/s^2', 'rad/s', 'rad/s'};
stribeckParMin = [0, 1, 0, 0, 0.02, 0.002];
stribeckParBaseMax = [30, 400, 80, 20, 20, 0.20];

%% ================================================================
% Setup
% ================================================================

diary(consoleLogFile);

fprintf('Data file:   %s\n', dataFile);
fprintf('Models path: %s\n', modelsFolder);
fprintf('Result file: %s\n\n', resultMatFile);

modelNamesToCheck = {viscousModelFile, coulombModelFile, stribeckModelFile};
for k = 1:numel(modelNamesToCheck)
    whichModel = which(modelNamesToCheck{k});
    if isempty(whichModel)
        error('MATLAB cannot find %s. Check modelsFolder and addpath(modelsFolder).', modelNamesToCheck{k});
    else
        fprintf('MATLAB will use %s from:\n  %s\n', modelNamesToCheck{k}, whichModel);
    end
end
fprintf('\n');

%% ================================================================
% Load data
% ================================================================

rawData = load(dataFile);
expData = struct();

for i = 1:nExp
    varName = sprintf('theta_2_%d_degrees', angles(i));
    data = rawData.(varName);

    expData(i).t = double(data.Time(:)); %#ok<SAGROW>
    expData(i).y_raw = double(data.Data(:));
    expData(i).u = zeros(size(expData(i).t));
    expData(i).initial_angle = angles(i);
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
        plot(expData(i).t, expData(i).y_raw, 'DisplayName', 'Raw data');
        hold on;
        xline(expData(i).t(startIndex(i)), '--', 'Start', 'HandleVisibility', 'off');
        xline(expData(i).t(endIndex(i)), '--', 'End', 'HandleVisibility', 'off');
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
        TsRaw = expData(i).t(2) - expData(i).t(1);
        fsRaw = 1 / TsRaw;

        [b, a] = butter(smoothingFilterOrder, smoothingCutoffHz/(fsRaw/2), 'low');
        expData(i).y_smooth = filtfilt(b, a, expData(i).y_raw);

        minDist = round(peakMinDistanceSeconds * fsRaw);
        [~, locMax] = findpeaks(expData(i).y_smooth, 'MinPeakDistance', minDist);
        [~, locMin] = findpeaks(-expData(i).y_smooth, 'MinPeakDistance', minDist);
        peakIndices = sort([locMax; locMin]);
        expData(i).peak_indices = peakIndices;

        nexttile;
        plot(expData(i).t, expData(i).y_raw, 'DisplayName', 'Raw data');
        hold on;
        plot(expData(i).t, expData(i).y_smooth, 'DisplayName', 'Smoothed data');
        plot(expData(i).t(peakIndices), expData(i).y_smooth(peakIndices), 'o', ...
            'DisplayName', 'Detected peaks');
        xline(expData(i).t(startIndex(i)), '--', 'Start', 'HandleVisibility', 'off');
        xline(expData(i).t(endIndex(i)), '--', 'End', 'HandleVisibility', 'off');

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
    expData(i).y_id = expData(i).y_raw(startIndex(i):endIndex(i));

    tailOffsetSamples = endIndex(i)+1:length(expData(i).y_raw);
    if isempty(tailOffsetSamples)
        offsetEstimate = mean(expData(i).y_raw(max(1, endIndex(i)-100):endIndex(i)));
    else
        offsetEstimate = mean(expData(i).y_raw(tailOffsetSamples));
    end

    expData(i).offsetEstimate = offsetEstimate;
    expData(i).y_id = expData(i).y_id - offsetEstimate;
    expData(i).t_id = expData(i).t(startIndex(i):endIndex(i));
    expData(i).t_id = expData(i).t_id - expData(i).t_id(1);
    expData(i).u_id = zeros(length(expData(i).y_id), 1);
end

Ts = expData(1).t_id(2) - expData(1).t_id(1);
fprintf('Sample time Ts = %.6f s\n\n', Ts);

%% ================================================================
% Plot sliced time-domain data
% ================================================================

if plotSlicedData
    figure('Name', 'Sliced and offset-corrected passive-link data');
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:nExp
        nexttile;
        plot(expData(i).t_id, expData(i).y_id);
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
        y = expData(i).y_id(:);
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
% Build tail-only data for low-speed Stribeck fitting
% ================================================================

for i = 1:nExp
    y = double(expData(i).y_id(:));
    t = double(expData(i).t_id(:));
    N = length(y);

    minTailSamples = max(20, round(minTailDuration / Ts));
    envelopeWindowSamples = max(5, round(envelopeWindowSeconds / Ts));

    env = movmax(abs(y), envelopeWindowSamples);
    candidateStart = find(env <= thetaTailMax, 1, 'first');

    if isempty(candidateStart)
        tailStart = round((1 - fallbackTailFraction) * N);
    else
        tailStart = candidateStart;
    end

    if N - tailStart + 1 < minTailSamples
        tailStart = max(1, N - minTailSamples + 1);
    end

    tailStart = max(1, min(tailStart, N));

    expData(i).tail_start_local_index = tailStart;
    expData(i).tail_start_time = t(tailStart);
    expData(i).y_tail = y(tailStart:end);
    expData(i).t_tail = t(tailStart:end) - t(tailStart);
    expData(i).u_tail = zeros(length(expData(i).y_tail), 1);
end

if plotTailSelection
    figure('Name', 'Tail sections for Stribeck identification');
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:nExp
        nexttile;
        plot(expData(i).t_id, expData(i).y_id, 'DisplayName', 'Full sliced data');
        hold on;
        plot(expData(i).t_tail + expData(i).tail_start_time, expData(i).y_tail, ...
            'LineWidth', 1.2, 'DisplayName', 'Tail used for Stribeck fit');
        xline(expData(i).tail_start_time, '--', 'Tail start', 'HandleVisibility', 'off');
        title(sprintf('%d degrees', angles(i)));
        xlabel('Time [s]');
        ylabel('\theta_2 [rad]');
        legend('Location', 'best');
    end

    sgtitle('Tail sections used for Stribeck identification');
    drawnow;
end

%% ================================================================
% Build iddata objects
% ================================================================

zFull = cell(nExp, 1);
zTail = cell(nExp, 1);

for i = 1:nExp
    zFull{i} = iddata(expData(i).y_id(:), expData(i).u_id(:), Ts);
    zFull{i}.Name = sprintf('%d_deg_full', angles(i));
    zFull{i}.InputName = {'u'};
    zFull{i}.InputUnit = {'-'};
    zFull{i}.OutputName = {'theta_2'};
    zFull{i}.OutputUnit = {'rad'};
    zFull{i}.TimeUnit = 's';

    zTail{i} = iddata(expData(i).y_tail(:), expData(i).u_tail(:), Ts);
    zTail{i}.Name = sprintf('%d_deg_tail', angles(i));
    zTail{i}.InputName = {'u'};
    zTail{i}.InputUnit = {'-'};
    zTail{i}.OutputName = {'theta_2'};
    zTail{i}.OutputUnit = {'rad'};
    zTail{i}.TimeUnit = 's';
end

zEst = zFull{idxEst(1)};
for k = 2:length(idxEst)
    zEst = merge(zEst, zFull{idxEst(k)});
end

zVal = zFull{idxVal(1)};
for k = 2:length(idxVal)
    zVal = merge(zVal, zFull{idxVal(k)});
end

zAll = zFull{idxAll(1)};
for k = 2:length(idxAll)
    zAll = merge(zAll, zFull{idxAll(k)});
end

zTailEst = zTail{idxEst(1)};
for k = 2:length(idxEst)
    zTailEst = merge(zTailEst, zTail{idxEst(k)});
end

zTailVal = zTail{idxVal(1)};
for k = 2:length(idxVal)
    zTailVal = merge(zTailVal, zTail{idxVal(k)});
end

zTailAll = zTail{idxAll(1)};
for k = 2:length(idxAll)
    zTailAll = merge(zTailAll, zTail{idxAll(k)});
end

%% ================================================================
% Estimate rough initial parameter guesses from the selected data
% ================================================================

p_b2_init = 0.067;
p_g2_init = 112;

try
    fPeaks = zeros(length(idxEst), 1);
    pbGuesses = zeros(length(idxEst), 1);

    for k = 1:length(idxEst)
        i = idxEst(k);
        y = double(expData(i).y_id(:));
        t = double(expData(i).t_id(:));

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
    p_g2_init = max(min(omega0^2, 400), 1);
    p_b2_init = max(min(median(pbGuesses), 30), 0.001);
catch ME
    warning('Could not compute data-derived initial guesses. Using fallback values instead: %s', ME.message);
end

fprintf('\nInitial rough guesses from data:\n');
fprintf('p_b2 = %.6f\n', p_b2_init);
fprintf('p_g2 = %.6f\n\n', p_g2_init);

%% ================================================================
% Estimation options
% ================================================================

compareOpt = compareOptions;
compareOpt.InitialCondition = 'estimate';

optSimple = nlgreyestOptions;
optSimple.Display = estimationDisplay;
optSimple.EstimateCovariance = estimateCovariance;
optSimple.SearchMethod = 'lm';
optSimple.SearchOptions.MaxIterations = maxIterationsSimple;

optCoulombMS = nlgreyestOptions;
optCoulombMS.Display = multistartDisplay;
optCoulombMS.EstimateCovariance = estimateCovariance;
optCoulombMS.SearchMethod = 'lm';
optCoulombMS.SearchOptions.MaxIterations = maxIterationsCoulombMultistart;

optStribeckMS = nlgreyestOptions;
optStribeckMS.Display = multistartDisplay;
optStribeckMS.EstimateCovariance = estimateCovariance;
optStribeckMS.SearchMethod = 'lm';
optStribeckMS.SearchOptions.MaxIterations = maxIterationsStribeckGrid;

%% ================================================================
% Stage 1: viscous + gravity model on full estimation data
% ================================================================

fprintf('\n============================================================\n');
fprintf('Stage 1: Viscous + gravity model, full data\n');
fprintf('============================================================\n');

% Initial states for estimation experiments, full-data segment.
theta0 = zeros(1, length(idxEst));
thetaDot0 = zeros(1, length(idxEst));
for k = 1:length(idxEst)
    i = idxEst(k);
    theta0(k) = double(expData(i).y_id(1));
    if estimateInitialSpeed
        tFitFull = double(expData(i).t_id(:));
        yFitFull = double(expData(i).y_id(:));
        nFit = min(velocityFitSamples, numel(tFitFull));
        if nFit < 2
            thetaDot0(k) = 0;
        else
            pFit = polyfit(tFitFull(1:nFit) - tFitFull(1), yFitFull(1:nFit), 1);
            thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1)));
        end
    else
        thetaDot0(k) = 0;
    end
end

mVisc0 = idnlgrey(viscousModelFile, [1 1 2], {p_b2_init; p_g2_init}, {theta0; thetaDot0}, 0);
mVisc0.Name = 'Passive rod viscous model';
mVisc0.TimeUnit = 's';
mVisc0.InputName = {'u'};
mVisc0.InputUnit = {'-'};
mVisc0.OutputName = {'theta_2'};
mVisc0.OutputUnit = {'rad'};
for k = 1:length(viscousParNames)
    mVisc0.Parameters(k).Name = viscousParNames{k};
    mVisc0.Parameters(k).Unit = viscousParUnits{k};
    mVisc0.Parameters(k).Minimum = viscousParMin(k);
    mVisc0.Parameters(k).Maximum = viscousParMax(k);
    mVisc0.Parameters(k).Fixed = false;
end
mVisc0.InitialStates(1).Name = 'theta_2';
mVisc0.InitialStates(1).Unit = 'rad';
mVisc0.InitialStates(2).Name = 'theta_2_dot';
mVisc0.InitialStates(2).Unit = 'rad/s';
mVisc0 = setinit(mVisc0, 'Fixed', {true(1, length(idxEst)); false(1, length(idxEst))});

mVisc = nlgreyest(zEst, mVisc0, optSimple);
mVisc.Name = 'Viscous';

fprintf('\nEstimated parameters for model: %s\n', mVisc.Name);
for k = 1:length(mVisc.Parameters)
    fprintf('%-12s = %12.6g', mVisc.Parameters(k).Name, mVisc.Parameters(k).Value);
    if mVisc.Parameters(k).Fixed
        fprintf('   fixed\n');
    else
        fprintf('   estimated\n');
    end
end

% Viscous validation and all-data models.
viscousParamValues = {mVisc.Parameters(1).Value; mVisc.Parameters(2).Value};

currentIdx = idxVal;
currentSegment = 'id';
theta0 = zeros(1, length(currentIdx));
thetaDot0 = zeros(1, length(currentIdx));
for k = 1:length(currentIdx)
    i = currentIdx(k);
    theta0(k) = double(expData(i).y_id(1));
    if estimateInitialSpeed
        tLocal = double(expData(i).t_id(:));
        yLocal = double(expData(i).y_id(:));
        nFit = min(velocityFitSamples, numel(tLocal));
        if nFit < 2
            thetaDot0(k) = 0;
        else
            pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1);
            thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1)));
        end
    else
        thetaDot0(k) = 0;
    end
end
mVisc_val = idnlgrey(viscousModelFile, [1 1 2], viscousParamValues, {theta0; thetaDot0}, 0);
mVisc_val.Name = mVisc.Name;
mVisc_val.TimeUnit = mVisc.TimeUnit;
mVisc_val.InputName = mVisc.InputName;
mVisc_val.InputUnit = mVisc.InputUnit;
mVisc_val.OutputName = mVisc.OutputName;
mVisc_val.OutputUnit = mVisc.OutputUnit;
for k = 1:length(mVisc.Parameters)
    mVisc_val.Parameters(k).Name = mVisc.Parameters(k).Name;
    mVisc_val.Parameters(k).Unit = mVisc.Parameters(k).Unit;
    mVisc_val.Parameters(k).Minimum = mVisc.Parameters(k).Minimum;
    mVisc_val.Parameters(k).Maximum = mVisc.Parameters(k).Maximum;
    mVisc_val.Parameters(k).Fixed = mVisc.Parameters(k).Fixed;
end
mVisc_val.InitialStates(1).Name = 'theta_2';
mVisc_val.InitialStates(1).Unit = 'rad';
mVisc_val.InitialStates(2).Name = 'theta_2_dot';
mVisc_val.InitialStates(2).Unit = 'rad/s';
mVisc_val = setinit(mVisc_val, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

currentIdx = idxAll;
theta0 = zeros(1, length(currentIdx));
thetaDot0 = zeros(1, length(currentIdx));
for k = 1:length(currentIdx)
    i = currentIdx(k);
    theta0(k) = double(expData(i).y_id(1));
    if estimateInitialSpeed
        tLocal = double(expData(i).t_id(:));
        yLocal = double(expData(i).y_id(:));
        nFit = min(velocityFitSamples, numel(tLocal));
        if nFit < 2
            thetaDot0(k) = 0;
        else
            pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1);
            thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1)));
        end
    else
        thetaDot0(k) = 0;
    end
end
mVisc_all = idnlgrey(viscousModelFile, [1 1 2], viscousParamValues, {theta0; thetaDot0}, 0);
mVisc_all.Name = mVisc.Name;
mVisc_all.TimeUnit = mVisc.TimeUnit;
mVisc_all.InputName = mVisc.InputName;
mVisc_all.InputUnit = mVisc.InputUnit;
mVisc_all.OutputName = mVisc.OutputName;
mVisc_all.OutputUnit = mVisc.OutputUnit;
for k = 1:length(mVisc.Parameters)
    mVisc_all.Parameters(k).Name = mVisc.Parameters(k).Name;
    mVisc_all.Parameters(k).Unit = mVisc.Parameters(k).Unit;
    mVisc_all.Parameters(k).Minimum = mVisc.Parameters(k).Minimum;
    mVisc_all.Parameters(k).Maximum = mVisc.Parameters(k).Maximum;
    mVisc_all.Parameters(k).Fixed = mVisc.Parameters(k).Fixed;
end
mVisc_all.InitialStates(1).Name = 'theta_2';
mVisc_all.InitialStates(1).Unit = 'rad';
mVisc_all.InitialStates(2).Name = 'theta_2_dot';
mVisc_all.InitialStates(2).Unit = 'rad/s';
mVisc_all = setinit(mVisc_all, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

[~, fitViscEst] = compare(zEst, mVisc, compareOpt);
[~, fitViscVal] = compare(zVal, mVisc_val, compareOpt);
[~, fitViscAll] = compare(zAll, mVisc_all, compareOpt);

fitVec = [];
if isnumeric(fitViscEst)
    fitVec = fitViscEst(:);
elseif iscell(fitViscEst)
    for k = 1:numel(fitViscEst)
        if isnumeric(fitViscEst{k}), fitVec = [fitVec; fitViscEst{k}(:)]; end %#ok<AGROW>
    end
end
fitVec = fitVec(isfinite(fitVec));
meanViscEst = mean(fitVec);

fitVec = [];
if isnumeric(fitViscVal)
    fitVec = fitViscVal(:);
elseif iscell(fitViscVal)
    for k = 1:numel(fitViscVal)
        if isnumeric(fitViscVal{k}), fitVec = [fitVec; fitViscVal{k}(:)]; end %#ok<AGROW>
    end
end
fitVec = fitVec(isfinite(fitVec));
meanViscVal = mean(fitVec);

fitVec = [];
if isnumeric(fitViscAll)
    fitVec = fitViscAll(:);
elseif iscell(fitViscAll)
    for k = 1:numel(fitViscAll)
        if isnumeric(fitViscAll{k}), fitVec = [fitVec; fitViscAll{k}(:)]; end %#ok<AGROW>
    end
end
fitVec = fitVec(isfinite(fitVec));
meanViscAll = mean(fitVec);

fprintf('\nViscous model mean estimation fit: %.2f %%\n', meanViscEst);
fprintf('Viscous model mean validation fit: %.2f %%\n', meanViscVal);
fprintf('Viscous model mean all-data fit: %.2f %%\n', meanViscAll);

%% ================================================================
% Stage 2: Coulomb-viscous model on full estimation data
% ================================================================

fprintf('\n============================================================\n');
fprintf('Stage 2: Coulomb-viscous model, full data\n');
fprintf('============================================================\n');

pCoul_b2 = mVisc.Parameters(1).Value;
pCoul_g2 = mVisc.Parameters(2).Value;
pCoul_c2 = 0.02 * pCoul_g2;
pCoul_eps = 0.03;

currentIdx = idxEst;
theta0 = zeros(1, length(currentIdx));
thetaDot0 = zeros(1, length(currentIdx));
for k = 1:length(currentIdx)
    i = currentIdx(k);
    theta0(k) = double(expData(i).y_id(1));
    if estimateInitialSpeed
        tLocal = double(expData(i).t_id(:));
        yLocal = double(expData(i).y_id(:));
        nFit = min(velocityFitSamples, numel(tLocal));
        if nFit < 2
            thetaDot0(k) = 0;
        else
            pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1);
            thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1)));
        end
    else
        thetaDot0(k) = 0;
    end
end

mCoul0 = idnlgrey(coulombModelFile, [1 1 2], {pCoul_b2; pCoul_g2; pCoul_c2; pCoul_eps}, {theta0; thetaDot0}, 0);
mCoul0.Name = 'Passive rod Coulomb-viscous model';
mCoul0.TimeUnit = 's';
mCoul0.InputName = {'u'};
mCoul0.InputUnit = {'-'};
mCoul0.OutputName = {'theta_2'};
mCoul0.OutputUnit = {'rad'};
for k = 1:length(coulombParNames)
    mCoul0.Parameters(k).Name = coulombParNames{k};
    mCoul0.Parameters(k).Unit = coulombParUnits{k};
    mCoul0.Parameters(k).Minimum = coulombParMin(k);
    mCoul0.Parameters(k).Maximum = coulombParMax(k);
    mCoul0.Parameters(k).Fixed = false;
end
mCoul0.Parameters(4).Fixed = true; % eps_v2 is fixed, as in the original workflow.
mCoul0.InitialStates(1).Name = 'theta_2';
mCoul0.InitialStates(1).Unit = 'rad';
mCoul0.InitialStates(2).Name = 'theta_2_dot';
mCoul0.InitialStates(2).Unit = 'rad/s';
mCoul0 = setinit(mCoul0, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

mCoulWarm = nlgreyest(zEst, mCoul0, optSimple);
mCoulWarm.Name = 'Coulomb-viscous warm start';

if runCoulombMultistart
    fprintf('\nRunning small Coulomb multistart on full estimation data...\n');
    bestCoulombScore = -Inf;
    mCoul = mCoulWarm;
    coulombResults = zeros(nStartsCoulomb, 4);
    rng(2);

    for sC = 1:nStartsCoulomb
        fprintf('Coulomb multistart %d / %d\n', sC, nStartsCoulomb);

        if sC == 1
            mCoulStart = mCoulWarm;
        else
            mCoulStart = mCoulWarm;
            for kp = 1:length(mCoulStart.Parameters)
                if ~mCoulStart.Parameters(kp).Fixed
                    current = mCoulStart.Parameters(kp).Value;
                    candidate = current * exp(0.45 * randn);
                    candidate = max(mCoulStart.Parameters(kp).Minimum, min(candidate, mCoulStart.Parameters(kp).Maximum));
                    mCoulStart.Parameters(kp).Value = candidate;
                end
            end
        end

        try
            mCoulTry = nlgreyest(zEst, mCoulStart, optCoulombMS);
            mCoulTry.Name = sprintf('Coulomb-viscous start %d', sC);

            coulombTryParamValues = cell(length(mCoulTry.Parameters), 1);
            for kp = 1:length(mCoulTry.Parameters)
                coulombTryParamValues{kp} = mCoulTry.Parameters(kp).Value;
            end

            currentIdx = idxVal;
            theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
            for kk = 1:length(currentIdx)
                ii = currentIdx(kk); theta0(kk) = expData(ii).y_id(1);
                tLocal = expData(ii).t_id(:); yLocal = expData(ii).y_id(:); nFit = min(velocityFitSamples, numel(tLocal));
                if estimateInitialSpeed && nFit >= 2
                    pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(kk) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1)));
                else
                    thetaDot0(kk) = 0;
                end
            end
            mCoulTryVal = idnlgrey(coulombModelFile, [1 1 2], coulombTryParamValues, {theta0; thetaDot0}, 0);
            mCoulTryVal.TimeUnit = mCoulTry.TimeUnit; mCoulTryVal.InputName = mCoulTry.InputName; mCoulTryVal.InputUnit = mCoulTry.InputUnit; mCoulTryVal.OutputName = mCoulTry.OutputName; mCoulTryVal.OutputUnit = mCoulTry.OutputUnit;
            for kp = 1:length(mCoulTry.Parameters)
                mCoulTryVal.Parameters(kp).Name = mCoulTry.Parameters(kp).Name; mCoulTryVal.Parameters(kp).Unit = mCoulTry.Parameters(kp).Unit; mCoulTryVal.Parameters(kp).Minimum = mCoulTry.Parameters(kp).Minimum; mCoulTryVal.Parameters(kp).Maximum = mCoulTry.Parameters(kp).Maximum; mCoulTryVal.Parameters(kp).Fixed = mCoulTry.Parameters(kp).Fixed;
            end
            mCoulTryVal.InitialStates(1).Name = 'theta_2'; mCoulTryVal.InitialStates(1).Unit = 'rad'; mCoulTryVal.InitialStates(2).Name = 'theta_2_dot'; mCoulTryVal.InitialStates(2).Unit = 'rad/s';
            mCoulTryVal = setinit(mCoulTryVal, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

            [~, fitTryEst] = compare(zEst, mCoulTry, compareOpt);
            [~, fitTryVal] = compare(zVal, mCoulTryVal, compareOpt);

            fitVec = [];
            if isnumeric(fitTryEst), fitVec = fitTryEst(:); elseif iscell(fitTryEst), for ff = 1:numel(fitTryEst), if isnumeric(fitTryEst{ff}), fitVec = [fitVec; fitTryEst{ff}(:)]; end, end, end %#ok<AGROW>
            fitVec = fitVec(isfinite(fitVec)); meanTryEst = mean(fitVec);
            fitVec = [];
            if isnumeric(fitTryVal), fitVec = fitTryVal(:); elseif iscell(fitTryVal), for ff = 1:numel(fitTryVal), if isnumeric(fitTryVal{ff}), fitVec = [fitVec; fitTryVal{ff}(:)]; end, end, end %#ok<AGROW>
            fitVec = fitVec(isfinite(fitVec)); meanTryVal = mean(fitVec);

            coulombScore = meanTryVal;
            coulombResults(sC, :) = [sC, meanTryEst, meanTryVal, NaN];
            fprintf('  Coulomb est: %.2f %% | val: %.2f %%\n', meanTryEst, meanTryVal);

            if coulombScore > bestCoulombScore
                bestCoulombScore = coulombScore;
                mCoul = mCoulTry;
            end
        catch ME
            warning('Coulomb multistart %d failed: %s', sC, ME.message);
            coulombResults(sC, :) = [sC, NaN, NaN, NaN];
        end
    end
    mCoul.Name = 'Best Coulomb-viscous model';
    coulombResultsTable = array2table(coulombResults, 'VariableNames', {'Start', 'MeanEstimationFit', 'MeanValidationFit', 'MeanAllDataFit'});
else
    mCoul = mCoulWarm;
    mCoul.Name = 'Coulomb-viscous';
    coulombResultsTable = table();
end

fprintf('\nEstimated parameters for model: %s\n', mCoul.Name);
for k = 1:length(mCoul.Parameters)
    fprintf('%-12s = %12.6g', mCoul.Parameters(k).Name, mCoul.Parameters(k).Value);
    if mCoul.Parameters(k).Fixed
        fprintf('   fixed\n');
    else
        fprintf('   estimated\n');
    end
end

% Coulomb validation and all-data models.
coulombParamValues = cell(length(mCoul.Parameters), 1);
for k = 1:length(mCoul.Parameters)
    coulombParamValues{k} = mCoul.Parameters(k).Value;
end

currentIdx = idxVal;
theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
for k = 1:length(currentIdx)
    i = currentIdx(k); theta0(k) = expData(i).y_id(1);
    if estimateInitialSpeed
        tLocal = expData(i).t_id(:); yLocal = expData(i).y_id(:); nFit = min(velocityFitSamples, numel(tLocal));
        if nFit < 2, thetaDot0(k) = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
    else
        thetaDot0(k) = 0;
    end
end
mCoul_val = idnlgrey(coulombModelFile, [1 1 2], coulombParamValues, {theta0; thetaDot0}, 0);
mCoul_val.Name = mCoul.Name; mCoul_val.TimeUnit = mCoul.TimeUnit; mCoul_val.InputName = mCoul.InputName; mCoul_val.InputUnit = mCoul.InputUnit; mCoul_val.OutputName = mCoul.OutputName; mCoul_val.OutputUnit = mCoul.OutputUnit;
for k = 1:length(mCoul.Parameters)
    mCoul_val.Parameters(k).Name = mCoul.Parameters(k).Name; mCoul_val.Parameters(k).Unit = mCoul.Parameters(k).Unit; mCoul_val.Parameters(k).Minimum = mCoul.Parameters(k).Minimum; mCoul_val.Parameters(k).Maximum = mCoul.Parameters(k).Maximum; mCoul_val.Parameters(k).Fixed = mCoul.Parameters(k).Fixed;
end
mCoul_val.InitialStates(1).Name = 'theta_2'; mCoul_val.InitialStates(1).Unit = 'rad'; mCoul_val.InitialStates(2).Name = 'theta_2_dot'; mCoul_val.InitialStates(2).Unit = 'rad/s';
mCoul_val = setinit(mCoul_val, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

currentIdx = idxAll;
theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
for k = 1:length(currentIdx)
    i = currentIdx(k); theta0(k) = expData(i).y_id(1);
    if estimateInitialSpeed
        tLocal = expData(i).t_id(:); yLocal = expData(i).y_id(:); nFit = min(velocityFitSamples, numel(tLocal));
        if nFit < 2, thetaDot0(k) = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
    else
        thetaDot0(k) = 0;
    end
end
mCoul_all = idnlgrey(coulombModelFile, [1 1 2], coulombParamValues, {theta0; thetaDot0}, 0);
mCoul_all.Name = mCoul.Name; mCoul_all.TimeUnit = mCoul.TimeUnit; mCoul_all.InputName = mCoul.InputName; mCoul_all.InputUnit = mCoul.InputUnit; mCoul_all.OutputName = mCoul.OutputName; mCoul_all.OutputUnit = mCoul.OutputUnit;
for k = 1:length(mCoul.Parameters)
    mCoul_all.Parameters(k).Name = mCoul.Parameters(k).Name; mCoul_all.Parameters(k).Unit = mCoul.Parameters(k).Unit; mCoul_all.Parameters(k).Minimum = mCoul.Parameters(k).Minimum; mCoul_all.Parameters(k).Maximum = mCoul.Parameters(k).Maximum; mCoul_all.Parameters(k).Fixed = mCoul.Parameters(k).Fixed;
end
mCoul_all.InitialStates(1).Name = 'theta_2'; mCoul_all.InitialStates(1).Unit = 'rad'; mCoul_all.InitialStates(2).Name = 'theta_2_dot'; mCoul_all.InitialStates(2).Unit = 'rad/s';
mCoul_all = setinit(mCoul_all, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

[~, fitCoulEst] = compare(zEst, mCoul, compareOpt);
[~, fitCoulVal] = compare(zVal, mCoul_val, compareOpt);
[~, fitCoulAll] = compare(zAll, mCoul_all, compareOpt);

fitVec = [];
if isnumeric(fitCoulEst), fitVec = fitCoulEst(:); elseif iscell(fitCoulEst), for k = 1:numel(fitCoulEst), if isnumeric(fitCoulEst{k}), fitVec = [fitVec; fitCoulEst{k}(:)]; end, end, end %#ok<AGROW>
fitVec = fitVec(isfinite(fitVec)); meanCoulEst = mean(fitVec);
fitVec = [];
if isnumeric(fitCoulVal), fitVec = fitCoulVal(:); elseif iscell(fitCoulVal), for k = 1:numel(fitCoulVal), if isnumeric(fitCoulVal{k}), fitVec = [fitVec; fitCoulVal{k}(:)]; end, end, end %#ok<AGROW>
fitVec = fitVec(isfinite(fitVec)); meanCoulVal = mean(fitVec);
fitVec = [];
if isnumeric(fitCoulAll), fitVec = fitCoulAll(:); elseif iscell(fitCoulAll), for k = 1:numel(fitCoulAll), if isnumeric(fitCoulAll{k}), fitVec = [fitVec; fitCoulAll{k}(:)]; end, end, end %#ok<AGROW>
fitVec = fitVec(isfinite(fitVec)); meanCoulAll = mean(fitVec);

fprintf('\nCoulomb model mean estimation fit: %.2f %%\n', meanCoulEst);
fprintf('Coulomb model mean validation fit: %.2f %%\n', meanCoulVal);
fprintf('Coulomb model mean all-data fit: %.2f %%\n', meanCoulAll);

%% ================================================================
% Stage 3: Grid search over fixed Stribeck velocity values
% ================================================================

fprintf('\n============================================================\n');
fprintf('Stage 3: Grid search over fixed Stribeck velocity values\n');
fprintf('============================================================\n');

bestModel = [];
bestTailValFit = -Inf;
bestSelectionScore = -Inf;

nGridRows = numel(stribeckVsCandidates) * nStartsStribeckPerVs;
results = zeros(nGridRows, 8);
row = 0;

rng(1);

for iVs = 1:numel(stribeckVsCandidates)
    vCandidate = stribeckVsCandidates(iVs);
    fprintf('\nTesting fixed v_s2 = %.4f rad/s\n', vCandidate);

    pGrid_b2 = mCoul.Parameters(1).Value;
    pGrid_g2 = mCoul.Parameters(2).Value;
    pGrid_c2 = mCoul.Parameters(3).Value;
    pGrid_sdelta = max(0.5 * pGrid_c2, 0.01);
    pGrid_vs = vCandidate;
    pGrid_eps = mCoul.Parameters(4).Value;

    currentIdx = idxEst;
    theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
    for k = 1:length(currentIdx)
        i = currentIdx(k); theta0(k) = expData(i).y_tail(1);
        if estimateInitialSpeed
            tLocal = expData(i).t_tail(:); yLocal = expData(i).y_tail(:); nFit = min(velocityFitSamples, numel(tLocal));
            if nFit < 2, thetaDot0(k) = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
        else
            thetaDot0(k) = 0;
        end
    end

    mGridBase = idnlgrey(stribeckModelFile, [1 1 2], {pGrid_b2; pGrid_g2; pGrid_c2; pGrid_sdelta; pGrid_vs; pGrid_eps}, {theta0; thetaDot0}, 0);
    mGridBase.Name = 'Passive rod report Stribeck model';
    mGridBase.TimeUnit = 's';
    mGridBase.InputName = {'u'};
    mGridBase.InputUnit = {'-'};
    mGridBase.OutputName = {'theta_2'};
    mGridBase.OutputUnit = {'rad'};
    stribeckParMax = stribeckParBaseMax;
    stribeckParMax(5) = max(2.50, 1.10 * pGrid_vs);
    for k = 1:length(stribeckParNames)
        mGridBase.Parameters(k).Name = stribeckParNames{k};
        mGridBase.Parameters(k).Unit = stribeckParUnits{k};
        mGridBase.Parameters(k).Minimum = stribeckParMin(k);
        mGridBase.Parameters(k).Maximum = stribeckParMax(k);
        mGridBase.Parameters(k).Fixed = false;
    end
    mGridBase.Parameters(1).Fixed = true;
    mGridBase.Parameters(2).Fixed = true;
    mGridBase.Parameters(3).Fixed = true;
    mGridBase.Parameters(4).Fixed = false;
    mGridBase.Parameters(5).Fixed = true;
    mGridBase.Parameters(6).Fixed = true;
    mGridBase.InitialStates(1).Name = 'theta_2'; mGridBase.InitialStates(1).Unit = 'rad';
    mGridBase.InitialStates(2).Name = 'theta_2_dot'; mGridBase.InitialStates(2).Unit = 'rad/s';
    mGridBase = setinit(mGridBase, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

    for sGrid = 1:nStartsStribeckPerVs
        row = row + 1;

        if sGrid == 1
            mStart = mGridBase;
        else
            mStart = mGridBase;
            current = mStart.Parameters(4).Value;
            candidate = current * exp(0.70 * randn);
            candidate = max(mStart.Parameters(4).Minimum, min(candidate, mStart.Parameters(4).Maximum));
            mStart.Parameters(4).Value = candidate;
        end

        try
            mTry = nlgreyest(zTailEst, mStart, optStribeckMS);
            mTry.Name = sprintf('Report Stribeck fixed vs %.4f start %d', vCandidate, sGrid);

            stribeckTryParams = cell(length(mTry.Parameters), 1);
            for kp = 1:length(mTry.Parameters)
                stribeckTryParams{kp} = mTry.Parameters(kp).Value;
            end

            % Tail validation model.
            currentIdx = idxVal;
            theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
            for kk = 1:length(currentIdx)
                ii = currentIdx(kk); theta0(kk) = expData(ii).y_tail(1);
                if estimateInitialSpeed
                    tLocal = expData(ii).t_tail(:); yLocal = expData(ii).y_tail(:); nFit = min(velocityFitSamples, numel(tLocal));
                    if nFit < 2, thetaDot0(kk) = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(kk) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
                else
                    thetaDot0(kk) = 0;
                end
            end
            mTry_tail_val = idnlgrey(stribeckModelFile, [1 1 2], stribeckTryParams, {theta0; thetaDot0}, 0);
            mTry_tail_val.Name = mTry.Name; mTry_tail_val.TimeUnit = mTry.TimeUnit; mTry_tail_val.InputName = mTry.InputName; mTry_tail_val.InputUnit = mTry.InputUnit; mTry_tail_val.OutputName = mTry.OutputName; mTry_tail_val.OutputUnit = mTry.OutputUnit;
            for kp = 1:length(mTry.Parameters)
                mTry_tail_val.Parameters(kp).Name = mTry.Parameters(kp).Name; mTry_tail_val.Parameters(kp).Unit = mTry.Parameters(kp).Unit; mTry_tail_val.Parameters(kp).Minimum = mTry.Parameters(kp).Minimum; mTry_tail_val.Parameters(kp).Maximum = mTry.Parameters(kp).Maximum; mTry_tail_val.Parameters(kp).Fixed = mTry.Parameters(kp).Fixed;
            end
            mTry_tail_val.InitialStates(1).Name = 'theta_2'; mTry_tail_val.InitialStates(1).Unit = 'rad'; mTry_tail_val.InitialStates(2).Name = 'theta_2_dot'; mTry_tail_val.InitialStates(2).Unit = 'rad/s';
            mTry_tail_val = setinit(mTry_tail_val, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

            % Full validation model.
            currentIdx = idxVal;
            theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
            for kk = 1:length(currentIdx)
                ii = currentIdx(kk); theta0(kk) = expData(ii).y_id(1);
                if estimateInitialSpeed
                    tLocal = expData(ii).t_id(:); yLocal = expData(ii).y_id(:); nFit = min(velocityFitSamples, numel(tLocal));
                    if nFit < 2, thetaDot0(kk) = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(kk) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
                else
                    thetaDot0(kk) = 0;
                end
            end
            mTry_full_val = idnlgrey(stribeckModelFile, [1 1 2], stribeckTryParams, {theta0; thetaDot0}, 0);
            mTry_full_val.Name = mTry.Name; mTry_full_val.TimeUnit = mTry.TimeUnit; mTry_full_val.InputName = mTry.InputName; mTry_full_val.InputUnit = mTry.InputUnit; mTry_full_val.OutputName = mTry.OutputName; mTry_full_val.OutputUnit = mTry.OutputUnit;
            for kp = 1:length(mTry.Parameters)
                mTry_full_val.Parameters(kp).Name = mTry.Parameters(kp).Name; mTry_full_val.Parameters(kp).Unit = mTry.Parameters(kp).Unit; mTry_full_val.Parameters(kp).Minimum = mTry.Parameters(kp).Minimum; mTry_full_val.Parameters(kp).Maximum = mTry.Parameters(kp).Maximum; mTry_full_val.Parameters(kp).Fixed = mTry.Parameters(kp).Fixed;
            end
            mTry_full_val.InitialStates(1).Name = 'theta_2'; mTry_full_val.InitialStates(1).Unit = 'rad'; mTry_full_val.InitialStates(2).Name = 'theta_2_dot'; mTry_full_val.InitialStates(2).Unit = 'rad/s';
            mTry_full_val = setinit(mTry_full_val, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

            % Full all-data model.
            currentIdx = idxAll;
            theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
            for kk = 1:length(currentIdx)
                ii = currentIdx(kk); theta0(kk) = expData(ii).y_id(1);
                if estimateInitialSpeed
                    tLocal = expData(ii).t_id(:); yLocal = expData(ii).y_id(:); nFit = min(velocityFitSamples, numel(tLocal));
                    if nFit < 2, thetaDot0(kk) = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(kk) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
                else
                    thetaDot0(kk) = 0;
                end
            end
            mTry_full_all = idnlgrey(stribeckModelFile, [1 1 2], stribeckTryParams, {theta0; thetaDot0}, 0);
            mTry_full_all.Name = mTry.Name; mTry_full_all.TimeUnit = mTry.TimeUnit; mTry_full_all.InputName = mTry.InputName; mTry_full_all.InputUnit = mTry.InputUnit; mTry_full_all.OutputName = mTry.OutputName; mTry_full_all.OutputUnit = mTry.OutputUnit;
            for kp = 1:length(mTry.Parameters)
                mTry_full_all.Parameters(kp).Name = mTry.Parameters(kp).Name; mTry_full_all.Parameters(kp).Unit = mTry.Parameters(kp).Unit; mTry_full_all.Parameters(kp).Minimum = mTry.Parameters(kp).Minimum; mTry_full_all.Parameters(kp).Maximum = mTry.Parameters(kp).Maximum; mTry_full_all.Parameters(kp).Fixed = mTry.Parameters(kp).Fixed;
            end
            mTry_full_all.InitialStates(1).Name = 'theta_2'; mTry_full_all.InitialStates(1).Unit = 'rad'; mTry_full_all.InitialStates(2).Name = 'theta_2_dot'; mTry_full_all.InitialStates(2).Unit = 'rad/s';
            mTry_full_all = setinit(mTry_full_all, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

            [~, fitTailEst] = compare(zTailEst, mTry, compareOpt);
            [~, fitTailVal] = compare(zTailVal, mTry_tail_val, compareOpt);
            [~, fitFullVal] = compare(zVal, mTry_full_val, compareOpt);
            [~, fitFullAll] = compare(zAll, mTry_full_all, compareOpt);

            fitVec = [];
            if isnumeric(fitTailEst), fitVec = fitTailEst(:); elseif iscell(fitTailEst), for ff = 1:numel(fitTailEst), if isnumeric(fitTailEst{ff}), fitVec = [fitVec; fitTailEst{ff}(:)]; end, end, end %#ok<AGROW>
            fitVec = fitVec(isfinite(fitVec)); meanTailEstFit = mean(fitVec);
            fitVec = [];
            if isnumeric(fitTailVal), fitVec = fitTailVal(:); elseif iscell(fitTailVal), for ff = 1:numel(fitTailVal), if isnumeric(fitTailVal{ff}), fitVec = [fitVec; fitTailVal{ff}(:)]; end, end, end %#ok<AGROW>
            fitVec = fitVec(isfinite(fitVec)); meanTailValFit = mean(fitVec);
            fitVec = [];
            if isnumeric(fitFullVal), fitVec = fitFullVal(:); elseif iscell(fitFullVal), for ff = 1:numel(fitFullVal), if isnumeric(fitFullVal{ff}), fitVec = [fitVec; fitFullVal{ff}(:)]; end, end, end %#ok<AGROW>
            fitVec = fitVec(isfinite(fitVec)); meanFullValFit = mean(fitVec);
            fitVec = [];
            if isnumeric(fitFullAll), fitVec = fitFullAll(:); elseif iscell(fitFullAll), for ff = 1:numel(fitFullAll), if isnumeric(fitFullAll{ff}), fitVec = [fitVec; fitFullAll{ff}(:)]; end, end, end %#ok<AGROW>
            fitVec = fitVec(isfinite(fitVec)); meanFullAllFit = mean(fitVec);

            selectionScore = tailValidationWeight * meanTailValFit + fullValidationWeight * meanFullValFit;
            pSdeltaValue = mTry.Parameters(4).Value;
            results(row, :) = [vCandidate, sGrid, pSdeltaValue, meanTailEstFit, meanTailValFit, meanFullValFit, meanFullAllFit, selectionScore];

            fprintf('  start %d/%d | p_sdelta2: %.4g | tail val: %.2f %% | full val: %.2f %% | score: %.2f\n', ...
                sGrid, nStartsStribeckPerVs, pSdeltaValue, meanTailValFit, meanFullValFit, selectionScore);

            if selectionScore > bestSelectionScore
                bestSelectionScore = selectionScore;
                bestTailValFit = meanTailValFit;
                bestModel = mTry;
            end
        catch ME
            warning('Stribeck grid failed for v_s2 %.4f, start %d: %s', vCandidate, sGrid, ME.message);
            results(row, :) = [vCandidate, sGrid, NaN, NaN, NaN, NaN, NaN, NaN];
        end
    end
end

if isempty(bestModel)
    error('All report/Stribeck grid runs failed. Check parameter bounds, tail data, and model file.');
else
    mStrBest = bestModel;
    mStrBest.Name = 'Best report Stribeck model, fixed-v_s2 grid';
end

fprintf('\nBest selected report/Stribeck model tail-validation fit: %.2f %%\n', bestTailValFit);
fprintf('\nEstimated parameters for model: %s\n', mStrBest.Name);
for k = 1:length(mStrBest.Parameters)
    fprintf('%-12s = %12.6g', mStrBest.Parameters(k).Name, mStrBest.Parameters(k).Value);
    if mStrBest.Parameters(k).Fixed
        fprintf('   fixed\n');
    else
        fprintf('   estimated\n');
    end
end

bestVs = mStrBest.Parameters(5).Value;
if abs(bestVs - max(stribeckVsCandidates)) < 1e-12
    warning(['Best Stribeck velocity v_s2 is at the upper grid boundary. ' ...
        'The optimum may lie above the tested grid. Extend stribeckVsCandidates ' ...
        'or treat this model as boundary-limited.']);
end

resultsTable = array2table(results, ...
    'VariableNames', {'v_s2', 'Start', 'p_sdelta2', 'MeanTailEstimationFit', ...
    'MeanTailValidationFit', 'MeanFullValidationFit', 'MeanFullAllDataFit', 'SelectionScore'});

disp(resultsTable);

%% ================================================================
% Final Stribeck evaluation
% ================================================================

stribeckBestParams = cell(length(mStrBest.Parameters), 1);
for k = 1:length(mStrBest.Parameters)
    stribeckBestParams{k} = mStrBest.Parameters(k).Value;
end

% Tail validation model.
currentIdx = idxVal;
theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
for k = 1:length(currentIdx)
    i = currentIdx(k); theta0(k) = expData(i).y_tail(1);
    if estimateInitialSpeed
        tLocal = expData(i).t_tail(:); yLocal = expData(i).y_tail(:); nFit = min(velocityFitSamples, numel(tLocal));
        if nFit < 2, thetaDot0(k) = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
    else
        thetaDot0(k) = 0;
    end
end
mStrBest_tail_val = idnlgrey(stribeckModelFile, [1 1 2], stribeckBestParams, {theta0; thetaDot0}, 0);
mStrBest_tail_val.Name = mStrBest.Name; mStrBest_tail_val.TimeUnit = mStrBest.TimeUnit; mStrBest_tail_val.InputName = mStrBest.InputName; mStrBest_tail_val.InputUnit = mStrBest.InputUnit; mStrBest_tail_val.OutputName = mStrBest.OutputName; mStrBest_tail_val.OutputUnit = mStrBest.OutputUnit;
for k = 1:length(mStrBest.Parameters)
    mStrBest_tail_val.Parameters(k).Name = mStrBest.Parameters(k).Name; mStrBest_tail_val.Parameters(k).Unit = mStrBest.Parameters(k).Unit; mStrBest_tail_val.Parameters(k).Minimum = mStrBest.Parameters(k).Minimum; mStrBest_tail_val.Parameters(k).Maximum = mStrBest.Parameters(k).Maximum; mStrBest_tail_val.Parameters(k).Fixed = mStrBest.Parameters(k).Fixed;
end
mStrBest_tail_val.InitialStates(1).Name = 'theta_2'; mStrBest_tail_val.InitialStates(1).Unit = 'rad'; mStrBest_tail_val.InitialStates(2).Name = 'theta_2_dot'; mStrBest_tail_val.InitialStates(2).Unit = 'rad/s';
mStrBest_tail_val = setinit(mStrBest_tail_val, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

% Full validation and all-data models use the last successful models from the grid style.
currentIdx = idxVal;
theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
for k = 1:length(currentIdx)
    i = currentIdx(k); theta0(k) = expData(i).y_id(1);
    if estimateInitialSpeed
        tLocal = expData(i).t_id(:); yLocal = expData(i).y_id(:); nFit = min(velocityFitSamples, numel(tLocal));
        if nFit < 2, thetaDot0(k) = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
    else
        thetaDot0(k) = 0;
    end
end
mStrBest_full_val = idnlgrey(stribeckModelFile, [1 1 2], stribeckBestParams, {theta0; thetaDot0}, 0);
mStrBest_full_val.Name = mStrBest.Name; mStrBest_full_val.TimeUnit = mStrBest.TimeUnit; mStrBest_full_val.InputName = mStrBest.InputName; mStrBest_full_val.InputUnit = mStrBest.InputUnit; mStrBest_full_val.OutputName = mStrBest.OutputName; mStrBest_full_val.OutputUnit = mStrBest.OutputUnit;
for k = 1:length(mStrBest.Parameters)
    mStrBest_full_val.Parameters(k).Name = mStrBest.Parameters(k).Name; mStrBest_full_val.Parameters(k).Unit = mStrBest.Parameters(k).Unit; mStrBest_full_val.Parameters(k).Minimum = mStrBest.Parameters(k).Minimum; mStrBest_full_val.Parameters(k).Maximum = mStrBest.Parameters(k).Maximum; mStrBest_full_val.Parameters(k).Fixed = mStrBest.Parameters(k).Fixed;
end
mStrBest_full_val.InitialStates(1).Name = 'theta_2'; mStrBest_full_val.InitialStates(1).Unit = 'rad'; mStrBest_full_val.InitialStates(2).Name = 'theta_2_dot'; mStrBest_full_val.InitialStates(2).Unit = 'rad/s';
mStrBest_full_val = setinit(mStrBest_full_val, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

currentIdx = idxAll;
theta0 = zeros(1, length(currentIdx)); thetaDot0 = zeros(1, length(currentIdx));
for k = 1:length(currentIdx)
    i = currentIdx(k); theta0(k) = expData(i).y_id(1);
    if estimateInitialSpeed
        tLocal = expData(i).t_id(:); yLocal = expData(i).y_id(:); nFit = min(velocityFitSamples, numel(tLocal));
        if nFit < 2, thetaDot0(k) = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0(k) = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
    else
        thetaDot0(k) = 0;
    end
end
mStrBest_full_all = idnlgrey(stribeckModelFile, [1 1 2], stribeckBestParams, {theta0; thetaDot0}, 0);
mStrBest_full_all.Name = mStrBest.Name; mStrBest_full_all.TimeUnit = mStrBest.TimeUnit; mStrBest_full_all.InputName = mStrBest.InputName; mStrBest_full_all.InputUnit = mStrBest.InputUnit; mStrBest_full_all.OutputName = mStrBest.OutputName; mStrBest_full_all.OutputUnit = mStrBest.OutputUnit;
for k = 1:length(mStrBest.Parameters)
    mStrBest_full_all.Parameters(k).Name = mStrBest.Parameters(k).Name; mStrBest_full_all.Parameters(k).Unit = mStrBest.Parameters(k).Unit; mStrBest_full_all.Parameters(k).Minimum = mStrBest.Parameters(k).Minimum; mStrBest_full_all.Parameters(k).Maximum = mStrBest.Parameters(k).Maximum; mStrBest_full_all.Parameters(k).Fixed = mStrBest.Parameters(k).Fixed;
end
mStrBest_full_all.InitialStates(1).Name = 'theta_2'; mStrBest_full_all.InitialStates(1).Unit = 'rad'; mStrBest_full_all.InitialStates(2).Name = 'theta_2_dot'; mStrBest_full_all.InitialStates(2).Unit = 'rad/s';
mStrBest_full_all = setinit(mStrBest_full_all, 'Fixed', {true(1, length(currentIdx)); false(1, length(currentIdx))});

[~, fitStrBestTailVal] = compare(zTailVal, mStrBest_tail_val, compareOpt);
[~, fitStrBestFullVal] = compare(zVal, mStrBest_full_val, compareOpt);
[~, fitStrBestFullAll] = compare(zAll, mStrBest_full_all, compareOpt);

fitVec = [];
if isnumeric(fitStrBestTailVal), fitVec = fitStrBestTailVal(:); elseif iscell(fitStrBestTailVal), for k = 1:numel(fitStrBestTailVal), if isnumeric(fitStrBestTailVal{k}), fitVec = [fitVec; fitStrBestTailVal{k}(:)]; end, end, end %#ok<AGROW>
fitVec = fitVec(isfinite(fitVec)); meanStrBestTailVal = mean(fitVec);
fitVec = [];
if isnumeric(fitStrBestFullVal), fitVec = fitStrBestFullVal(:); elseif iscell(fitStrBestFullVal), for k = 1:numel(fitStrBestFullVal), if isnumeric(fitStrBestFullVal{k}), fitVec = [fitVec; fitStrBestFullVal{k}(:)]; end, end, end %#ok<AGROW>
fitVec = fitVec(isfinite(fitVec)); meanStrBestFullVal = mean(fitVec);
fitVec = [];
if isnumeric(fitStrBestFullAll), fitVec = fitStrBestFullAll(:); elseif iscell(fitStrBestFullAll), for k = 1:numel(fitStrBestFullAll), if isnumeric(fitStrBestFullAll{k}), fitVec = [fitVec; fitStrBestFullAll{k}(:)]; end, end, end %#ok<AGROW>
fitVec = fitVec(isfinite(fitVec)); meanStrBestFullAll = mean(fitVec);

fprintf('\nBest report/Stribeck model mean tail validation fit: %.2f %%\n', meanStrBestTailVal);
fprintf('Best report/Stribeck model mean full validation fit: %.2f %%\n', meanStrBestFullVal);
fprintf('Best report/Stribeck model mean full all-data fit: %.2f %%\n', meanStrBestFullAll);

%% ================================================================
% Parameter uncertainty and identifiability diagnostics
% ================================================================

modelList = {mVisc, mCoul, mStrBest};
modelLabels = {'Viscous', 'Coulomb-viscous', 'Best Stribeck fixed-vs grid'};
allParameterUncertainty = table();
identifiabilitySummary = table();
parameterCorrelation = struct();

for im = 1:numel(modelList)
    mDiag = modelList{im};
    modelLabel = modelLabels{im};

    nPar = length(mDiag.Parameters);
    model = repmat({modelLabel}, nPar, 1);
    parameter = cell(nPar, 1);
    value = zeros(nPar, 1);
    fixed = false(nPar, 1);
    variance = NaN(nPar, 1);
    stdDev = NaN(nPar, 1);
    relativeStdPercent = NaN(nPar, 1);
    ci95Lower = NaN(nPar, 1);
    ci95Upper = NaN(nPar, 1);

    for kp = 1:nPar
        parameter{kp} = mDiag.Parameters(kp).Name;
        value(kp) = double(mDiag.Parameters(kp).Value);
        fixed(kp) = mDiag.Parameters(kp).Fixed;
    end

    freeIdx = find(~fixed);
    try
        covFree = getcov(mDiag, 'value', 'free');
        if ~isempty(covFree) && isnumeric(covFree)
            if size(covFree, 1) == length(freeIdx)
                variance(freeIdx) = diag(covFree);
            elseif size(covFree, 1) == nPar
                variance = diag(covFree);
            end
        end
    catch ME
        warning('Could not extract free-parameter covariance for model %s: %s', modelLabel, ME.message);
        try
            covAll = getcov(mDiag, 'value');
            if ~isempty(covAll) && isnumeric(covAll) && size(covAll, 1) == nPar
                variance = diag(covAll);
                variance(fixed) = NaN;
            end
        catch ME2
            warning('Could not extract full parameter covariance for model %s: %s', modelLabel, ME2.message);
        end
    end

    validVariance = isfinite(variance);
    variance(validVariance) = max(variance(validVariance), 0);
    stdDev(validVariance) = sqrt(variance(validVariance));

    for kp = 1:nPar
        if isfinite(stdDev(kp))
            ci95Lower(kp) = value(kp) - 1.96 * stdDev(kp);
            ci95Upper(kp) = value(kp) + 1.96 * stdDev(kp);
            if abs(value(kp)) > eps
                relativeStdPercent(kp) = 100 * stdDev(kp) / abs(value(kp));
            end
        end
    end

    tbl = table(model, parameter, value, fixed, variance, stdDev, relativeStdPercent, ci95Lower, ci95Upper, ...
        'VariableNames', {'Model', 'Parameter', 'Value', 'Fixed', 'Variance', 'StdDev', 'RelativeStdPercent', 'CI95Lower', 'CI95Upper'});
    allParameterUncertainty = [allParameterUncertainty; tbl]; %#ok<AGROW>

    freeRows = ~tbl.Fixed;
    numParameters = height(tbl);
    numFreeParameters = sum(freeRows);
    if any(freeRows)
        maxRelativeStdPercent = max(tbl.RelativeStdPercent(freeRows), [], 'omitnan');
    else
        maxRelativeStdPercent = NaN;
    end

    maxAbsCorrelation = NaN;
    covarianceConditionNumber = NaN;
    corrTbl = table();
    try
        covFree = getcov(mDiag, 'value', 'free');
        if ~isempty(covFree) && isnumeric(covFree) && size(covFree, 1) >= 2
            stdFree = sqrt(max(diag(covFree), 0));
            denominator = stdFree * stdFree.';
            corrFree = covFree ./ denominator;
            corrFree(denominator == 0) = NaN;
            corrFree(1:size(corrFree,1)+1:end) = 1;
            mask = ~eye(size(corrFree));
            maxAbsCorrelation = max(abs(corrFree(mask)), [], 'omitnan');
            covarianceConditionNumber = cond(covFree);

            freeNames = parameter(freeIdx);
            rowNames = matlab.lang.makeUniqueStrings(freeNames);
            varNames = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(freeNames));
            corrTbl = array2table(corrFree, 'VariableNames', varNames, 'RowNames', rowNames);
        end
    catch
        maxAbsCorrelation = NaN;
        covarianceConditionNumber = NaN;
        corrTbl = table();
    end

    if im == 1
        parameterCorrelation.viscous = corrTbl;
    elseif im == 2
        parameterCorrelation.coulomb = corrTbl;
    else
        parameterCorrelation.stribeckBest = corrTbl;
    end

    summaryRow = table({modelLabel}, numParameters, numFreeParameters, maxRelativeStdPercent, maxAbsCorrelation, covarianceConditionNumber, ...
        'VariableNames', {'Model', 'NumParameters', 'NumFreeParameters', 'MaxRelativeStdPercent', 'MaxAbsFreeParameterCorrelation', 'FreeCovarianceConditionNumber'});
    identifiabilitySummary = [identifiabilitySummary; summaryRow]; %#ok<AGROW>
end

fprintf('\nParameter uncertainty table:\n');
disp(allParameterUncertainty);

fprintf('\nIdentifiability summary:\n');
disp(identifiabilitySummary);

fprintf('\nBest Stribeck free-parameter correlation matrix:\n');
disp(parameterCorrelation.stribeckBest);

writetable(allParameterUncertainty, uncertaintyCsvFile);
writetable(identifiabilitySummary, identifiabilityCsvFile);

%% ================================================================
% Final model-vs-measurement plots
% ================================================================

if plotFinalFits
    plotModels = {mVisc, mCoul, mStrBest};
    plotModelFiles = {viscousModelFile, coulombModelFile, stribeckModelFile};
    plotTitles = {'Viscous model, full-data fit', ...
                  'Coulomb-viscous model, full-data fit', ...
                  'Best report Stribeck model, fixed-v_s2 grid, shown on full data'};

    for im = 1:numel(plotModels)
        mPlotSource = plotModels{im};
        modelFilePlot = plotModelFiles{im};
        figure('Name', plotTitles{im});
        tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        for k = 1:nExp
            theta0Single = expData(k).y_id(1);
            if estimateInitialSpeed
                tLocal = expData(k).t_id(:); yLocal = expData(k).y_id(:); nFit = min(velocityFitSamples, numel(tLocal));
                if nFit < 2, thetaDot0Single = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0Single = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
            else
                thetaDot0Single = 0;
            end

            paramsSingle = cell(length(mPlotSource.Parameters), 1);
            for kp = 1:length(mPlotSource.Parameters)
                paramsSingle{kp} = mPlotSource.Parameters(kp).Value;
            end
            mPlot = idnlgrey(modelFilePlot, [1 1 2], paramsSingle, {theta0Single; thetaDot0Single}, 0);
            mPlot.Name = mPlotSource.Name; mPlot.TimeUnit = mPlotSource.TimeUnit; mPlot.InputName = mPlotSource.InputName; mPlot.InputUnit = mPlotSource.InputUnit; mPlot.OutputName = mPlotSource.OutputName; mPlot.OutputUnit = mPlotSource.OutputUnit;
            for kp = 1:length(mPlotSource.Parameters)
                mPlot.Parameters(kp).Name = mPlotSource.Parameters(kp).Name; mPlot.Parameters(kp).Unit = mPlotSource.Parameters(kp).Unit; mPlot.Parameters(kp).Minimum = mPlotSource.Parameters(kp).Minimum; mPlot.Parameters(kp).Maximum = mPlotSource.Parameters(kp).Maximum; mPlot.Parameters(kp).Fixed = mPlotSource.Parameters(kp).Fixed;
            end
            mPlot.InitialStates(1).Name = 'theta_2'; mPlot.InitialStates(1).Unit = 'rad'; mPlot.InitialStates(2).Name = 'theta_2_dot'; mPlot.InitialStates(2).Unit = 'rad/s';
            mPlot = setinit(mPlot, 'Fixed', {true; false});

            [yHatPlot, fitPlot] = compare(zFull{k}, mPlot, compareOpt);
            yRaw = get(yHatPlot, 'OutputData');
            if iscell(yRaw), yRaw = yRaw{1}; end
            yModel = double(yRaw(:));
            measured = expData(k).y_id(:);
            n = min(length(measured), length(yModel));

            nexttile;
            plot(expData(k).t_id(1:n), measured(1:n), 'DisplayName', 'Measured');
            hold on;
            plot(expData(k).t_id(1:n), yModel(1:n), '--', 'DisplayName', 'Model');
            if isnumeric(fitPlot) && isfinite(fitPlot(1))
                fitText = sprintf(', fit %.1f%%', fitPlot(1));
            else
                fitText = '';
            end
            title(sprintf('%d degrees%s', angles(k), fitText));
            xlabel('Time [s]');
            ylabel('\theta_2 [rad]');
            legend('Location', 'best');
        end
        sgtitle(plotTitles{im});
    end

    % Tail-only view for checking the low-speed Stribeck region directly.
    figure('Name', 'Best Stribeck model, tail-only view');
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:nExp
        theta0Single = expData(k).y_tail(1);
        if estimateInitialSpeed
            tLocal = expData(k).t_tail(:); yLocal = expData(k).y_tail(:); nFit = min(velocityFitSamples, numel(tLocal));
            if nFit < 2, thetaDot0Single = 0; else, pFit = polyfit(tLocal(1:nFit) - tLocal(1), yLocal(1:nFit), 1); thetaDot0Single = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, pFit(1))); end
        else
            thetaDot0Single = 0;
        end

        paramsSingle = cell(length(mStrBest.Parameters), 1);
        for kp = 1:length(mStrBest.Parameters)
            paramsSingle{kp} = mStrBest.Parameters(kp).Value;
        end
        mPlot = idnlgrey(stribeckModelFile, [1 1 2], paramsSingle, {theta0Single; thetaDot0Single}, 0);
        mPlot.Name = mStrBest.Name; mPlot.TimeUnit = mStrBest.TimeUnit; mPlot.InputName = mStrBest.InputName; mPlot.InputUnit = mStrBest.InputUnit; mPlot.OutputName = mStrBest.OutputName; mPlot.OutputUnit = mStrBest.OutputUnit;
        for kp = 1:length(mStrBest.Parameters)
            mPlot.Parameters(kp).Name = mStrBest.Parameters(kp).Name; mPlot.Parameters(kp).Unit = mStrBest.Parameters(kp).Unit; mPlot.Parameters(kp).Minimum = mStrBest.Parameters(kp).Minimum; mPlot.Parameters(kp).Maximum = mStrBest.Parameters(kp).Maximum; mPlot.Parameters(kp).Fixed = mStrBest.Parameters(kp).Fixed;
        end
        mPlot.InitialStates(1).Name = 'theta_2'; mPlot.InitialStates(1).Unit = 'rad'; mPlot.InitialStates(2).Name = 'theta_2_dot'; mPlot.InitialStates(2).Unit = 'rad/s';
        mPlot = setinit(mPlot, 'Fixed', {true; false});

        [yHatPlot, fitPlot] = compare(zTail{k}, mPlot, compareOpt);
        yRaw = get(yHatPlot, 'OutputData');
        if iscell(yRaw), yRaw = yRaw{1}; end
        yModel = double(yRaw(:));
        measured = expData(k).y_tail(:);
        n = min(length(measured), length(yModel));

        nexttile;
        plot(expData(k).t_tail(1:n), measured(1:n), 'DisplayName', 'Measured tail');
        hold on;
        plot(expData(k).t_tail(1:n), yModel(1:n), '--', 'DisplayName', 'Model');
        if isnumeric(fitPlot) && isfinite(fitPlot(1))
            fitText = sprintf(', fit %.1f%%', fitPlot(1));
        else
            fitText = '';
        end
        title(sprintf('%d degrees%s', angles(k), fitText));
        xlabel('Tail time [s]');
        ylabel('\theta_2 [rad]');
        legend('Location', 'best');
    end
    sgtitle('Best report Stribeck model, fixed-v_s2 grid, tail-only view');
    drawnow;
end

%% ================================================================
% Save results
% ================================================================

passiveID.mVisc = mVisc;
passiveID.mCoul = mCoul;
passiveID.mStrBest = mStrBest;
passiveID.coulombResultsTable = coulombResultsTable;
passiveID.resultsTable = resultsTable;
passiveID.idxEst = idxEst;
passiveID.idxVal = idxVal;
passiveID.tailSettings.thetaMax = thetaTailMax;
passiveID.tailSettings.minTailDuration = minTailDuration;
passiveID.tailSettings.fallbackFraction = fallbackTailFraction;
passiveID.tailSettings.envelopeWindow = envelopeWindowSeconds;
passiveID.stribeckVsCandidates = stribeckVsCandidates;
passiveID.nStartsStribeckPerVs = nStartsStribeckPerVs;
passiveID.stribeckSelection.tailValidationWeight = tailValidationWeight;
passiveID.stribeckSelection.fullValidationWeight = fullValidationWeight;
passiveID.parameterUncertainty = allParameterUncertainty;
passiveID.identifiabilitySummary = identifiabilitySummary;
passiveID.parameterCorrelation = parameterCorrelation;
passiveID.expData = expData;
passiveID.config.dataFile = dataFile;
passiveID.config.modelsFolder = modelsFolder;
passiveID.config.figurePosition = figurePosition;
passiveID.config.startIndex = startIndex;
passiveID.config.endIndex = endIndex;

save(resultMatFile, 'passiveID');

fprintf('\nSaved result file: %s\n', resultMatFile);
fprintf('Saved uncertainty table: %s\n', uncertaintyCsvFile);
fprintf('Saved identifiability table: %s\n', identifiabilityCsvFile);

diary off;
