clear;          % Removes all variables from the workspace
clear functions;  % Clears cached local/helper functions so MATLAB cannot use an older version
clc;            % Clears the Command Window
close all;      % Closes all figure windows
fprintf('Running rod_2_identification_v11_estimated_initial_speed.m\n');

%% Load data
raw_data = load("identification_data\rod_2_identification\matlab.mat");

%% Unpack the loaded data into a struct
angles = 15:15:120;

for i = 1:length(angles)
    var_name = sprintf('theta_2_%d_degrees', angles(i));
    data = raw_data.(var_name);

    expData(i).t = double(data.Time(:));
    expData(i).y_raw = double(data.Data(:));
    expData(i).u = zeros(size(expData(i).t));
    expData(i).initial_angle = angles(i);
end

%% Plot raw data in time domain
figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(angles)
    nexttile;
    plot(expData(i).t, expData(i).y_raw, "DisplayName", "Raw Data");
    title(sprintf('%d Degrees', angles(i)));
    ylabel("Angle (rad)");
    xlabel("Time (s)");
end

%% Plot data with peak detection to help select useful truncation indices
figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(angles)

    Ts_raw = expData(i).t(2) - expData(i).t(1);
    fs_raw = 1/Ts_raw;

    % Filter only for peak detection. Identification uses unfiltered data.
    fc = 5;
    [b, a] = butter(3, fc/(fs_raw/2), 'low');

    expData(i).y_smooth = filtfilt(b, a, expData(i).y_raw);

    minDist = round(0.3 * fs_raw);
    [~, locMax] = findpeaks(expData(i).y_smooth, 'MinPeakDistance', minDist);
    [~, locMin] = findpeaks(-expData(i).y_smooth, 'MinPeakDistance', minDist);
    peak_indices = sort([locMax; locMin]);

    nexttile;
    plot(expData(i).t, expData(i).y_raw, 'DisplayName', 'Raw data');
    hold on;
    plot(expData(i).t, expData(i).y_smooth, 'DisplayName', 'Smoothed data');
    plot(expData(i).t(peak_indices), expData(i).y_smooth(peak_indices), 'o', ...
        'DisplayName', 'Detected peaks');

    legend;
    title(sprintf('Data %d Degrees', angles(i)));
    xlabel('Time (s)');
    ylabel('Angle (rad)');
end

%% Set indices for truncated data
expData(1).start_index = 919;
expData(1).end_index = 1703;

expData(2).start_index = 1078;
expData(2).end_index = 2594;

expData(3).start_index = 1089;
expData(3).end_index = 3339;

expData(4).start_index = 955;
expData(4).end_index = 3643;

expData(5).start_index = 1091;
expData(5).end_index = 3969;

expData(6).start_index = 1051;
expData(6).end_index = 4656;

expData(7).start_index = 851;
expData(7).end_index = 4857;

expData(8).start_index = 802;
expData(8).end_index = 4984;

for i = 1:length(angles)
    expData(i).y_id = expData(i).y_raw(expData(i).start_index:expData(i).end_index);

    tailOffsetSamples = expData(i).end_index + 1:length(expData(i).y_raw);

    if isempty(tailOffsetSamples)
        offsetEstimate = mean(expData(i).y_raw(max(1, expData(i).end_index - 100):expData(i).end_index));
    else
        offsetEstimate = mean(expData(i).y_raw(tailOffsetSamples));
    end

    expData(i).y_id = expData(i).y_id - offsetEstimate;

    expData(i).t_id = expData(i).t(expData(i).start_index:expData(i).end_index);
    expData(i).t_id = expData(i).t_id - expData(i).t_id(1);

    expData(i).u_id = zeros(length(expData(i).y_id), 1);
end

%% Plot truncated data in time domain
figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(angles)
    nexttile;
    plot(expData(i).t_id, expData(i).y_id);
    title(sprintf('Truncated Raw Data %d Degrees', angles(i)));
    xlabel('Time (s)');
    ylabel('Angle (rad)');
end

%% Plot truncated data in frequency domain
figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(angles)

    Ts_fft = expData(i).t_id(2) - expData(i).t_id(1);
    fs_fft = 1/Ts_fft;

    data_length = length(expData(i).y_id);

    y_fft = fft(expData(i).y_id);
    P2 = abs(y_fft / data_length);
    P1 = P2(1:floor(data_length/2)+1);
    P1(2:end-1) = 2 * P1(2:end-1);

    f = fs_fft * (0:floor(data_length/2)) / data_length;

    nexttile;
    plot(f, P1, 'LineWidth', 1.2);

    title(sprintf('Raw Data %d Degrees', angles(i)));
    xlabel("Frequency (Hz)");
    ylabel("Amplitude");
    xscale('log');
end

%% ================================================================
% Identification setup - passive rod 2, report-derived Stribeck model
%
% This version keeps the simple full-data models for reference and then
% tests the friction model derived in the report:
%
%   F_Str(theta_dot) = p_b2 theta_dot
%       + [p_c2 + p_sdelta exp(-(theta_dot/v_s)^2)] tanh(theta_dot/eps_v)
%
% where:
%   p_c2     = kinetic Coulomb friction strength
%   p_sdelta = static-minus-Coulomb friction strength, constrained >= 0
%   v_s      = Stribeck velocity
%
% Because the measured data is shifted so the final resting angle is zero,
% the passive-link gravitational term is implemented as
%
%   -p_g2 sin(theta)
%
% without a separate output/angle offset parameter.
%
% Staged identification:
%   1) Viscous + gravity model on full selected data
%   2) Coulomb-viscous model on full selected data
%   3) Grid search over fixed Stribeck velocities; only p_sdelta is estimated
%
% In all stages theta_dot(0) is estimated as an initial condition to reduce
% phase/lead-lag errors from imperfect manual truncation.
% ================================================================

%% Settings
idxEst = [1 3 5 7];  % 15, 45, 75, 105 degrees
idxVal = [2 4 6 8];  % 30, 60, 90, 120 degrees
idxAll = 1:length(angles);

% Tail selection focuses the Stribeck fit on the small-angle, low-speed region.
tailSettings.thetaMax = 0.30;          % rad, low-amplitude region
tailSettings.minTailDuration = 3.0;    % s
tailSettings.fallbackFraction = 0.40;  % last 40% if amplitude logic gives too little data
tailSettings.envelopeWindow = 0.60;    % s moving maximum window

% Runtime control.
% The Coulomb model is simple and the output showed that multistart barely
% improves it. Use one estimation from the viscous result to save runtime.
runCoulombMultistart = false;
nStartsCoulomb = 1;

% The Stribeck parameters p_sdelta and v_s are strongly correlated if both
% are freely estimated from free-decay tail data. Therefore v_s is treated as
% a grid-searched hyperparameter, and only p_sdelta is estimated for each v_s.
% The previous optimum kept moving to the upper grid boundary. Therefore the
% grid is extended further upward. If the selected value still lands at the
% largest candidate, the result should be treated as boundary-limited and the
% grid can be extended again.
stribeckVsCandidates = [2.50 2.75 3 3.25 3.50];

% Only p_sdelta is free for a fixed v_s, so repeated random starts give almost
% identical results. One start per candidate is enough and much faster.
nStartsStribeckPerVs = 1;

% Final Stribeck model selection. Full-validation fit is weighted more heavily
% to avoid choosing a tail-only solution that stops the full oscillation too early.
stribeckSelection.tailValidationWeight = 0.30;
stribeckSelection.fullValidationWeight = 0.70;

% Initial angular velocity handling.
% The truncation indices are close to peaks, but not necessarily exactly at
% zero angular velocity. If theta_dot(0) is forced to zero while the data
% starts slightly before/after a peak, the model can appear to lead or lag the
% measurement even when the friction model is good. Therefore theta(0) is kept
% fixed to the measured first sample, while theta_dot(0) is estimated as a
% nuisance initial condition during nlgreyest. The initial guess for theta_dot(0)
% is obtained from a short linear fit to the first samples of each segment.
initialStateSettings.estimateInitialSpeed = true;
initialStateSettings.velocityFitSamples = 9;
initialStateSettings.maxAbsInitialSpeed = 5.0;  % rad/s, only used for the initial guess clamp

% During compare/plotting, estimate the initial condition again. This makes the
% validation plots test the model dynamics without penalising a small manually
% selected start-index timing error. If this is not available in your MATLAB
% version, the helper automatically falls back to normal compare(...).
useEstimatedInitialConditionForCompare = true; %#ok<NASGU>

% The model files are written automatically so this script is self-contained.
overwritePassiveModelFiles = true;
scriptFolder = getCurrentScriptFolder();
writePassiveModelFiles(scriptFolder, overwritePassiveModelFiles);
addpath(scriptFolder);

%% Build tail-only data for low-speed friction fitting
expData = buildTailData(expData, angles, tailSettings);

figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(angles)
    nexttile;
    plot(expData(i).t_id, expData(i).y_id, 'DisplayName', 'Full identification data');
    hold on;
    plot(expData(i).t_tail + expData(i).tail_start_time, expData(i).y_tail, ...
        'LineWidth', 1.2, 'DisplayName', 'Tail used for Stribeck fit');
    xline(expData(i).tail_start_time, '--', 'Tail start', 'HandleVisibility', 'off');
    title(sprintf('%d Degrees', angles(i)));
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('Location', 'best');
end

sgtitle('Tail sections used for report/Stribeck identification');

%% Create iddata objects
Ts = expData(1).t_id(2) - expData(1).t_id(1);

zFull = createIddataCell(expData, angles, Ts, 'id');
zTail = createIddataCell(expData, angles, Ts, 'tail');

zEst = mergeIddataExperiments(zFull, idxEst);
zVal = mergeIddataExperiments(zFull, idxVal);
zAll = mergeIddataExperiments(zFull, idxAll);

zTailEst = mergeIddataExperiments(zTail, idxEst);
zTailVal = mergeIddataExperiments(zTail, idxVal);
zTailAll = mergeIddataExperiments(zTail, idxAll);

%% Estimate rough initial guesses from the data
[pb2_init, pg2_init] = estimateInitialPassiveGuess(expData, idxEst);

fprintf("\nInitial rough guesses from data:\n");
fprintf("p_b2 = %.6f\n", pb2_init);
fprintf("p_g2 = %.6f\n\n", pg2_init);

%% Estimation options
optSimple = nlgreyestOptions;
optSimple.Display = 'on';
optSimple.EstimateCovariance = true;
optSimple.SearchOptions.MaxIterations = 60;

optStribeckWarm = nlgreyestOptions;
optStribeckWarm.Display = 'on';
optStribeckWarm.EstimateCovariance = true;
optStribeckWarm.SearchOptions.MaxIterations = 80;

optCoulombMS = nlgreyestOptions;
optCoulombMS.Display = 'off';
optCoulombMS.EstimateCovariance = true;
optCoulombMS.SearchOptions.MaxIterations = 50;

optStribeckMS = nlgreyestOptions;
optStribeckMS.Display = 'off';
optStribeckMS.EstimateCovariance = true;
optStribeckMS.SearchOptions.MaxIterations = 90;

%% ================================================================
% Stage 1: viscous + gravity model on full estimation data
% ================================================================

fprintf("\n============================================================\n");
fprintf("Stage 1: Viscous + gravity model, full data\n");
fprintf("============================================================\n");

pVisc0.p_b2 = pb2_init;
pVisc0.p_g2 = pg2_init;

mVisc0 = createPassiveViscousModel(expData, idxEst, pVisc0, 'id');
mVisc = nlgreyest(zEst, mVisc0, optSimple);
mVisc.Name = 'Viscous';

dispPassiveParameters(mVisc);
uncertaintyVisc = parameterUncertaintyTable(mVisc, 'Viscous');
disp(uncertaintyVisc);

mVisc_val = setPassiveInitialStates(mVisc, expData, idxVal, 'id');

[~, fitViscEst] = comparePassiveModel(zEst, mVisc);
[~, fitViscVal] = comparePassiveModel(zVal, mVisc_val);
[~, fitViscAll] = comparePassiveModel(zAll, setPassiveInitialStates(mVisc, expData, idxAll, 'id'));

fprintf("\nViscous model mean estimation fit: %.2f %%\n", meanFitValue(fitViscEst));
fprintf("Viscous model mean validation fit: %.2f %%\n", meanFitValue(fitViscVal));
fprintf("Viscous model mean all-data fit: %.2f %%\n", meanFitValue(fitViscAll));

%% ================================================================
% Stage 2: Coulomb-viscous model on full estimation data
% ================================================================

fprintf("\n============================================================\n");
fprintf("Stage 2: Coulomb-viscous model, full data\n");
fprintf("============================================================\n");

pCoul0.p_b2  = getParValue(mVisc, 'p_b2');
pCoul0.p_g2  = getParValue(mVisc, 'p_g2');
pCoul0.p_c2  = 0.02 * pCoul0.p_g2;
pCoul0.eps_v = 0.03;

mCoul0 = createPassiveCoulombModel(expData, idxEst, pCoul0, 'id');

% eps_v mostly controls the smooth sign approximation around zero speed.
% Fixing it here improves convergence speed.
mCoul0 = setParFixed(mCoul0, 'eps_v', true);

mCoulWarm = nlgreyest(zEst, mCoul0, optSimple);
mCoulWarm.Name = 'Coulomb-viscous warm start';

if runCoulombMultistart
    fprintf("\nRunning small Coulomb multistart on full estimation data...\n");

    mCoulBase = mCoulWarm;
    mCoulBase = setParFixed(mCoulBase, 'p_b2', false);
    mCoulBase = setParFixed(mCoulBase, 'p_g2', false);
    mCoulBase = setParFixed(mCoulBase, 'p_c2', false);
    mCoulBase = setParFixed(mCoulBase, 'eps_v', true);

    bestCoulombScore = -Inf;
    mCoul = mCoulWarm;
    coulombResults = zeros(nStartsCoulomb, 4);

    rng(2);

    for sC = 1:nStartsCoulomb
        fprintf("Coulomb multistart %d / %d\n", sC, nStartsCoulomb);

        if sC == 1
            mCoulStart = mCoulBase;
        else
            mCoulStart = randomizeFreeParameters(mCoulBase, 0.45);
        end

        try
            mCoulTry = nlgreyest(zEst, mCoulStart, optCoulombMS);
            mCoulTry.Name = sprintf('Coulomb-viscous start %d', sC);

            mCoulTry_val = setPassiveInitialStates(mCoulTry, expData, idxVal, 'id');
            mCoulTry_all = setPassiveInitialStates(mCoulTry, expData, idxAll, 'id');

            [~, fitCoulTryEst] = comparePassiveModel(zEst, mCoulTry);
            [~, fitCoulTryVal] = comparePassiveModel(zVal, mCoulTry_val);
            [~, fitCoulTryAll] = comparePassiveModel(zAll, mCoulTry_all);

            meanCoulTryEst = meanFitValue(fitCoulTryEst);
            meanCoulTryVal = meanFitValue(fitCoulTryVal);
            meanCoulTryAll = meanFitValue(fitCoulTryAll);

            % Select the Coulomb model using validation first, with all-data
            % fit as a mild tie-breaker. This avoids an overfitted baseline
            % friction estimate being inherited by the Stribeck model.
            coulombScore = meanCoulTryVal + 0.10 * meanCoulTryAll;

            coulombResults(sC, :) = [sC, meanCoulTryEst, meanCoulTryVal, meanCoulTryAll];

            fprintf("  Coulomb est: %.2f %% | val: %.2f %% | all: %.2f %%\n", ...
                meanCoulTryEst, meanCoulTryVal, meanCoulTryAll);

            if coulombScore > bestCoulombScore
                bestCoulombScore = coulombScore;
                mCoul = mCoulTry;
            end

        catch ME
            warning("Coulomb multistart %d failed: %s", sC, ME.message);
            coulombResults(sC, :) = [sC, NaN, NaN, NaN];
        end
    end

    mCoul.Name = 'Best Coulomb-viscous model';
    coulombResultsTable = array2table(coulombResults, ...
        'VariableNames', {'Start', 'MeanEstimationFit', 'MeanValidationFit', 'MeanAllDataFit'});
    disp(coulombResultsTable);
else
    mCoul = mCoulWarm;
    mCoul.Name = 'Coulomb-viscous';
    coulombResultsTable = table();
end

dispPassiveParameters(mCoul);
uncertaintyCoul = parameterUncertaintyTable(mCoul, 'Coulomb-viscous');
disp(uncertaintyCoul);

mCoul_val = setPassiveInitialStates(mCoul, expData, idxVal, 'id');

[~, fitCoulEst] = comparePassiveModel(zEst, mCoul);
[~, fitCoulVal] = comparePassiveModel(zVal, mCoul_val);
[~, fitCoulAll] = comparePassiveModel(zAll, setPassiveInitialStates(mCoul, expData, idxAll, 'id'));

fprintf("\nCoulomb model mean estimation fit: %.2f %%\n", meanFitValue(fitCoulEst));
fprintf("Coulomb model mean validation fit: %.2f %%\n", meanFitValue(fitCoulVal));
fprintf("Coulomb model mean all-data fit: %.2f %%\n", meanFitValue(fitCoulAll));

%% ================================================================
% Stage 3: Grid search over fixed Stribeck velocity values
%
% For each candidate v_s, only p_sdelta is estimated. This produces a much
% better-conditioned identification problem than estimating p_sdelta and v_s
% simultaneously from free-decay tail data.
% ================================================================

fprintf("\n============================================================\n");
fprintf("Stage 3: Grid search over fixed Stribeck velocity values\n");
fprintf("============================================================\n");

bestModel = [];
bestTailValFit = -Inf;
bestSelectionScore = -Inf;

nGridRows = numel(stribeckVsCandidates) * nStartsStribeckPerVs;
results = zeros(nGridRows, 8);  % [v_s, start, p_sdelta, tail_est, tail_val, full_val, full_all, score]
row = 0;

rng(1);

for iVs = 1:numel(stribeckVsCandidates)

    vCandidate = stribeckVsCandidates(iVs);

    fprintf("\nTesting fixed v_s = %.4f rad/s\n", vCandidate);

    pGrid.p_b2 = getParValue(mCoul, 'p_b2');
    pGrid.p_g2 = getParValue(mCoul, 'p_g2');
    pGrid.p_c2 = getParValue(mCoul, 'p_c2');
    pGrid.p_sdelta = max(0.5 * pGrid.p_c2, 0.01);
    pGrid.v_s = vCandidate;
    pGrid.eps_v = getParValue(mCoul, 'eps_v');

    mGridBase = createPassiveStribeckModel(expData, idxEst, pGrid, 'tail');
    mGridBase = setParFixed(mGridBase, 'p_b2', true);
    mGridBase = setParFixed(mGridBase, 'p_g2', true);
    mGridBase = setParFixed(mGridBase, 'p_c2', true);
    mGridBase = setParFixed(mGridBase, 'eps_v', true);
    mGridBase = setParFixed(mGridBase, 'v_s', true);
    mGridBase = setParFixed(mGridBase, 'p_sdelta', false);

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

            mTry_tail_val = setPassiveInitialStates(mTry, expData, idxVal, 'tail');
            mTry_full_val = setPassiveInitialStates(mTry, expData, idxVal, 'id');
            mTry_all_full = setPassiveInitialStates(mTry, expData, idxAll, 'id');

            [~, fitTailEst] = comparePassiveModel(zTailEst, mTry);
            [~, fitTailVal] = comparePassiveModel(zTailVal, mTry_tail_val);
            [~, fitFullVal] = comparePassiveModel(zVal, mTry_full_val);
            [~, fitFullAll] = comparePassiveModel(zAll, mTry_all_full);

            meanTailEstFit = meanFitValue(fitTailEst);
            meanTailValFit = meanFitValue(fitTailVal);
            meanFullValFit = meanFitValue(fitFullVal);
            meanFullAllFit = meanFitValue(fitFullAll);

            selectionScore = ...
                stribeckSelection.tailValidationWeight * meanTailValFit + ...
                stribeckSelection.fullValidationWeight * meanFullValFit;

            pSdeltaValue = getParValue(mTry, 'p_sdelta');
            results(row, :) = [vCandidate, sGrid, pSdeltaValue, meanTailEstFit, ...
                meanTailValFit, meanFullValFit, meanFullAllFit, selectionScore];

            fprintf("  start %d/%d | p_sdelta: %.4g | tail val: %.2f %% | full val: %.2f %% | score: %.2f\n", ...
                sGrid, nStartsStribeckPerVs, pSdeltaValue, meanTailValFit, meanFullValFit, selectionScore);

            if selectionScore > bestSelectionScore
                bestSelectionScore = selectionScore;
                bestTailValFit = meanTailValFit;
                bestModel = mTry;
            end

        catch ME
            warning("Stribeck grid failed for v_s %.4f, start %d: %s", vCandidate, sGrid, ME.message);
            results(row, :) = [vCandidate, sGrid, NaN, NaN, NaN, NaN, NaN, NaN];
        end
    end
end

if isempty(bestModel)
    error("All report/Stribeck grid runs failed. Check parameter bounds, tail data, and model file.");
else
    mStrBest = bestModel;
    mStrBest.Name = 'Best report Stribeck model, fixed-v_s grid';
end

fprintf("\nBest selected report/Stribeck model tail-validation fit: %.2f %%\n", bestTailValFit);
dispPassiveParameters(mStrBest);
uncertaintyStrBest = parameterUncertaintyTable(mStrBest, 'Best Stribeck fixed-vs grid');
disp(uncertaintyStrBest);

bestVs = getParValue(mStrBest, 'v_s');
if abs(bestVs - max(stribeckVsCandidates)) < 1e-12
    warning(['Best Stribeck velocity v_s is at the upper grid boundary. ' ...
        'The optimum may lie above the tested grid. Extend stribeckVsCandidates ' ...
        'or treat this model as boundary-limited.']);
end

resultsTable = array2table(results, ...
    'VariableNames', {'v_s', 'Start', 'p_sdelta', 'MeanTailEstimationFit', ...
    'MeanTailValidationFit', 'MeanFullValidationFit', 'MeanFullAllDataFit', 'SelectionScore'});

disp(resultsTable);

%% ================================================================
% Final model evaluation and plots for all eight measurements
% ================================================================

mStrBest_tail_val = setPassiveInitialStates(mStrBest, expData, idxVal, 'tail');
mStrBest_full_val = setPassiveInitialStates(mStrBest, expData, idxVal, 'id');
mStrBest_full_all = setPassiveInitialStates(mStrBest, expData, idxAll, 'id');

[~, fitStrBestTailVal] = comparePassiveModel(zTailVal, mStrBest_tail_val);
[~, fitStrBestFullVal] = comparePassiveModel(zVal, mStrBest_full_val);
[~, fitStrBestFullAll] = comparePassiveModel(zAll, mStrBest_full_all);

fprintf("\nBest report/Stribeck model mean tail validation fit: %.2f %%\n", meanFitValue(fitStrBestTailVal));
fprintf("Best report/Stribeck model mean full validation fit: %.2f %%\n", meanFitValue(fitStrBestFullVal));
fprintf("Best report/Stribeck model mean full all-data fit: %.2f %%\n", meanFitValue(fitStrBestFullAll));

%% ================================================================
% Parameter uncertainty and identifiability diagnostics
% ================================================================

allParameterUncertainty = [
    uncertaintyVisc;
    uncertaintyCoul;
    uncertaintyStrBest
];

identifiabilitySummary = [
    parameterIdentifiabilitySummary(mVisc, 'Viscous');
    parameterIdentifiabilitySummary(mCoul, 'Coulomb-viscous');
    parameterIdentifiabilitySummary(mStrBest, 'Best Stribeck fixed-vs grid')
];

correlationVisc = parameterCorrelationTable(mVisc);
correlationCoul = parameterCorrelationTable(mCoul);
correlationStrBest = parameterCorrelationTable(mStrBest);

fprintf("\nParameter uncertainty table:\n");
disp(allParameterUncertainty);

fprintf("\nIdentifiability summary:\n");
disp(identifiabilitySummary);

fprintf("\nBest Stribeck free-parameter correlation matrix:\n");
disp(correlationStrBest);

writetable(allParameterUncertainty, 'parameter_uncertainty_summary.csv');
writetable(identifiabilitySummary, 'parameter_identifiability_summary.csv');

% One 4x2 figure per model. Each figure shows all eight experiments.
plotModelAgainstAllMeasurements(mVisc, 'Viscous model, full-data fit', zFull, expData, angles, 'id', true);
plotModelAgainstAllMeasurements(mCoul, 'Coulomb-viscous model, full-data fit', zFull, expData, angles, 'id', true);
plotModelAgainstAllMeasurements(mStrBest, 'Best report Stribeck model, fixed-v_s grid, shown on full data', zFull, expData, angles, 'id', true);

% Tail-only view for checking the low-speed region directly.
plotModelAgainstAllMeasurements(mStrBest, 'Best report Stribeck model, fixed-v_s grid, tail-only view', zTail, expData, angles, 'tail', false);

%% Save results
passiveID.mVisc = mVisc;
passiveID.mCoul = mCoul;
passiveID.mStrBest = mStrBest;
passiveID.coulombResultsTable = coulombResultsTable;
passiveID.resultsTable = resultsTable;
passiveID.idxEst = idxEst;
passiveID.idxVal = idxVal;
passiveID.tailSettings = tailSettings;
passiveID.stribeckVsCandidates = stribeckVsCandidates;
passiveID.nStartsStribeckPerVs = nStartsStribeckPerVs;
passiveID.stribeckSelection = stribeckSelection;
passiveID.parameterUncertainty = allParameterUncertainty;
passiveID.identifiabilitySummary = identifiabilitySummary;
passiveID.parameterCorrelation.viscous = correlationVisc;
passiveID.parameterCorrelation.coulomb = correlationCoul;
passiveID.parameterCorrelation.stribeckBest = correlationStrBest;
passiveID.expData = expData;

save('passive_rod_2_identification_report_stribeck_result.mat', 'passiveID');

%% ================================================================
% Local helper functions
% ================================================================

function folder = getCurrentScriptFolder()
%GETCURRENTSCRIPTFOLDER Return folder of this script, or pwd if unavailable.

fullName = mfilename('fullpath');

if isempty(fullName)
    folder = pwd;
else
    folder = fileparts(fullName);
end

if isempty(folder)
    folder = pwd;
end

end

function expData = buildTailData(expData, angles, tailSettings)
%BUILDTailDATA Build small-amplitude tail sections for low-speed fitting.

for i = 1:length(angles)

    y = double(expData(i).y_id(:));
    t = double(expData(i).t_id(:));
    Ts = t(2) - t(1);
    N = length(y);

    minTailSamples = max(20, round(tailSettings.minTailDuration / Ts));
    envelopeWindowSamples = max(5, round(tailSettings.envelopeWindow / Ts));

    env = movmax(abs(y), envelopeWindowSamples);
    candidateStart = find(env <= tailSettings.thetaMax, 1, 'first');

    if isempty(candidateStart)
        tailStart = round((1 - tailSettings.fallbackFraction) * N);
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

end

function z = createIddataCell(expData, angles, Ts, segment)
%CREATEIDDATACELL Create one iddata object per experiment.

z = cell(length(angles), 1);

for i = 1:length(angles)

    switch segment
        case 'id'
            y = double(expData(i).y_id(:));
            u = double(expData(i).u_id(:));
        case 'tail'
            y = double(expData(i).y_tail(:));
            u = double(expData(i).u_tail(:));
        otherwise
            error("Unknown segment '%s'. Use 'id' or 'tail'.", segment);
    end

    if length(y) ~= length(u)
        error("Experiment %d, segment %s: input and output lengths do not match.", i, segment);
    end

    z{i} = iddata(y, u, Ts);
    z{i}.Name = sprintf('%d_deg_%s', angles(i), segment);
    z{i}.InputName = {'u'};
    z{i}.OutputName = {'theta'};
    z{i}.InputUnit = {'-'};
    z{i}.OutputUnit = {'rad'};
    z{i}.TimeUnit = 's';
end

end

function zMerged = mergeIddataExperiments(z, idx)
%MERGEIDDATAEXPERIMENTS Merge selected iddata experiments.

zMerged = z{idx(1)};

for k = 2:length(idx)
    zMerged = merge(zMerged, z{idx(k)});
end

end

function [pb2_init, pg2_init] = estimateInitialPassiveGuess(expData, idxEst)
%ESTIMATEINITIALPASSIVEGUESS Estimate rough p_b2 and p_g2 guesses.

f_peaks = zeros(length(idxEst), 1);
pb_guesses = zeros(length(idxEst), 1);

for k = 1:length(idxEst)

    i = idxEst(k);

    y = double(expData(i).y_id(:));
    t = double(expData(i).t_id(:));

    Ts = t(2) - t(1);
    fs = 1 / Ts;
    N = length(y);

    y0 = y - mean(y);

    Y = fft(y0);
    P2 = abs(Y / N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2 * P1(2:end-1);

    f = fs * (0:floor(N/2)) / N;

    valid = f > 0.2 & f < 10;

    if ~any(valid)
        error("No valid frequency range found for experiment %d.", i);
    end

    [~, idxMax] = max(P1(valid));
    fValid = f(valid);
    f_peaks(k) = fValid(idxMax);

    omega_d = 2 * pi * f_peaks(k);

    minPeakDist = round(0.35 / Ts);
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

    pb_guesses(k) = 2 * zeta * omega_d;
end

omega0 = 2 * pi * median(f_peaks);
pg2_init = omega0^2;

pb2_init = median(pb_guesses);

pg2_init = max(min(pg2_init, 400), 1);
pb2_init = max(min(pb2_init, 30), 0.001);

end

function m = createPassiveViscousModel(expData, idx, p0, segment)
%CREATEPASSIVEVISCOUSMODEL Create idnlgrey viscous passive-rod model.

FileName = 'passive_rod_viscous_model';
Order = [1 1 2];
TsModel = 0;

Parameters = {p0.p_b2; p0.p_g2};
InitialStates = passiveInitialStateCell(expData, idx, segment);

m = idnlgrey(FileName, Order, Parameters, InitialStates, TsModel);

m.Name = 'Passive rod viscous model';
m.TimeUnit = 's';

m.InputName = {'u'};
m.InputUnit = {'-'};
m.OutputName = {'theta'};
m.OutputUnit = {'rad'};

m.Parameters(1).Name = 'p_b2';
m.Parameters(1).Unit = '1/s';
m.Parameters(1).Minimum = 0;
m.Parameters(1).Maximum = 30;
m.Parameters(1).Fixed = false;

m.Parameters(2).Name = 'p_g2';
m.Parameters(2).Unit = '1/s^2';
m.Parameters(2).Minimum = 1;
m.Parameters(2).Maximum = 400;
m.Parameters(2).Fixed = false;

m = fixPassiveInitialStates(m);

end

function m = createPassiveCoulombModel(expData, idx, p0, segment)
%CREATEPASSIVECOULOMBMODEL Create idnlgrey Coulomb-viscous passive-rod model.

FileName = 'passive_rod_coulomb_model';
Order = [1 1 2];
TsModel = 0;

Parameters = {p0.p_b2; p0.p_g2; p0.p_c2; p0.eps_v};
InitialStates = passiveInitialStateCell(expData, idx, segment);

m = idnlgrey(FileName, Order, Parameters, InitialStates, TsModel);

m.Name = 'Passive rod Coulomb-viscous model';
m.TimeUnit = 's';

m.InputName = {'u'};
m.InputUnit = {'-'};
m.OutputName = {'theta'};
m.OutputUnit = {'rad'};

m.Parameters(1).Name = 'p_b2';
m.Parameters(1).Unit = '1/s';
m.Parameters(1).Minimum = 0;
m.Parameters(1).Maximum = 30;
m.Parameters(1).Fixed = false;

m.Parameters(2).Name = 'p_g2';
m.Parameters(2).Unit = '1/s^2';
m.Parameters(2).Minimum = 1;
m.Parameters(2).Maximum = 400;
m.Parameters(2).Fixed = false;

m.Parameters(3).Name = 'p_c2';
m.Parameters(3).Unit = 'rad/s^2';
m.Parameters(3).Minimum = 0;
m.Parameters(3).Maximum = 80;
m.Parameters(3).Fixed = false;

m.Parameters(4).Name = 'eps_v';
m.Parameters(4).Unit = 'rad/s';
m.Parameters(4).Minimum = 0.002;
m.Parameters(4).Maximum = 0.20;
m.Parameters(4).Fixed = false;

m = fixPassiveInitialStates(m);

end

function m = createPassiveStribeckModel(expData, idx, p0, segment)
%CREATEPASSIVESTRIBECKMODEL Create report-derived Stribeck passive-rod model.

FileName = 'passive_rod_stribeck_report_model';
Order = [1 1 2];
TsModel = 0;

Parameters = {p0.p_b2; p0.p_g2; p0.p_c2; p0.p_sdelta; p0.v_s; p0.eps_v};
InitialStates = passiveInitialStateCell(expData, idx, segment);

m = idnlgrey(FileName, Order, Parameters, InitialStates, TsModel);

m.Name = 'Passive rod report Stribeck model';
m.TimeUnit = 's';

m.InputName = {'u'};
m.InputUnit = {'-'};
m.OutputName = {'theta'};
m.OutputUnit = {'rad'};

m.Parameters(1).Name = 'p_b2';
m.Parameters(1).Unit = '1/s';
m.Parameters(1).Minimum = 0;
m.Parameters(1).Maximum = 30;
m.Parameters(1).Fixed = false;

m.Parameters(2).Name = 'p_g2';
m.Parameters(2).Unit = '1/s^2';
m.Parameters(2).Minimum = 1;
m.Parameters(2).Maximum = 400;
m.Parameters(2).Fixed = false;

m.Parameters(3).Name = 'p_c2';
m.Parameters(3).Unit = 'rad/s^2';
m.Parameters(3).Minimum = 0;
m.Parameters(3).Maximum = 80;
m.Parameters(3).Fixed = false;

% p_sdelta = F_S - F_C, so F_S >= F_C is enforced by p_sdelta >= 0.
m.Parameters(4).Name = 'p_sdelta';
m.Parameters(4).Unit = 'rad/s^2';
m.Parameters(4).Minimum = 0;
m.Parameters(4).Maximum = 20;
m.Parameters(4).Fixed = false;

m.Parameters(5).Name = 'v_s';
m.Parameters(5).Unit = 'rad/s';
m.Parameters(5).Minimum = 0.02;
% Make the upper bound depend on the candidate value. This avoids errors
% when the grid is extended above the previous hard-coded maximum.
m.Parameters(5).Maximum = max(2.50, 1.10 * p0.v_s);
m.Parameters(5).Fixed = false;

m.Parameters(6).Name = 'eps_v';
m.Parameters(6).Unit = 'rad/s';
m.Parameters(6).Minimum = 0.002;
m.Parameters(6).Maximum = 0.20;
m.Parameters(6).Fixed = false;

m = fixPassiveInitialStates(m);

end

function X0 = passiveInitialStateCell(expData, idx, segment)
%PASSIVEINITIALSTATECELL Create initial states for multi-experiment idnlgrey.
%
% theta(0) is fixed to the measured first sample. theta_dot(0) is initialised
% from the first few samples of the selected segment and then estimated as a
% nuisance initial condition by nlgreyest. This avoids forcing the experiment
% to start exactly at a turning point.

nExp = length(idx);

theta0 = zeros(1, nExp);
thetaDot0 = zeros(1, nExp);

for k = 1:nExp

    i = idx(k);

    switch segment
        case 'id'
            theta0(k) = double(expData(i).y_id(1));
        case 'tail'
            theta0(k) = double(expData(i).y_tail(1));
        otherwise
            error("Unknown segment '%s'. Use 'id' or 'tail'.", segment);
    end

    thetaDot0(k) = estimateInitialAngularVelocity(expData, i, segment);
end

X0 = {theta0; thetaDot0};

end

function thetaDot0 = estimateInitialAngularVelocity(expData, i, segment)
%ESTIMATEINITIALANGULARVELOCITY Estimate an initial angular velocity guess.
%
% The estimate is only an initial guess. Since theta_dot(0) is free during
% nlgreyest, the optimiser can adjust it. A short linear fit is less noisy than
% a single finite difference.

switch segment
    case 'id'
        t = double(expData(i).t_id(:));
        y = double(expData(i).y_id(:));
    case 'tail'
        t = double(expData(i).t_tail(:));
        y = double(expData(i).y_tail(:));
    otherwise
        error("Unknown segment '%s'. Use 'id' or 'tail'.", segment);
end

if numel(t) < 2 || numel(y) < 2
    thetaDot0 = 0;
    return;
end

nFit = min(9, numel(t));
tFit = t(1:nFit) - t(1);
yFit = y(1:nFit);

if all(abs(tFit) < eps)
    thetaDot0 = 0;
else
    p = polyfit(tFit, yFit, 1);
    thetaDot0 = p(1);
end

% Clamp only the starting guess to avoid extreme values caused by noisy or very
% short windows. The parameter itself is still estimated freely by nlgreyest.
maxAbsInitialSpeed = 5.0;
thetaDot0 = max(-maxAbsInitialSpeed, min(maxAbsInitialSpeed, thetaDot0));

end

function mOut = setPassiveInitialStates(mIn, expData, idx, segment)
%SETPASSIVEINITIALSTATES Return a model with initial states for idx/segment.
%
% This function deliberately rebuilds the idnlgrey model instead of trying to
% resize InitialStates on the existing model. MATLAB does not allow resizing
% those InitialStates arrays after the object has been created.

mOut = rebuildPassiveModelWithInitialStates(mIn, expData, idx, segment);

end

function mOut = rebuildPassiveModelWithInitialStates(mIn, expData, idx, segment)
%REBUILDPASSIVEMODELWITHINITIALSTATES Recreate passive model for new experiment count.

nPar = length(mIn.Parameters);

switch nPar
    case 2
        p0.p_b2 = getParValue(mIn, 'p_b2');
        p0.p_g2 = getParValue(mIn, 'p_g2');
        mOut = createPassiveViscousModel(expData, idx, p0, segment);

    case 4
        p0.p_b2  = getParValue(mIn, 'p_b2');
        p0.p_g2  = getParValue(mIn, 'p_g2');
        p0.p_c2  = getParValue(mIn, 'p_c2');
        p0.eps_v = getParValue(mIn, 'eps_v');
        mOut = createPassiveCoulombModel(expData, idx, p0, segment);

    case 6
        p0.p_b2 = getParValue(mIn, 'p_b2');
        p0.p_g2 = getParValue(mIn, 'p_g2');
        p0.p_c2 = getParValue(mIn, 'p_c2');
        p0.p_sdelta = getParValue(mIn, 'p_sdelta');
        p0.v_s = getParValue(mIn, 'v_s');
        p0.eps_v = getParValue(mIn, 'eps_v');
        mOut = createPassiveStribeckModel(expData, idx, p0, segment);

    otherwise
        error('Unsupported passive model with %d parameters.', nPar);
end

% Preserve parameter settings from the estimated model.
for k = 1:nPar
    mOut.Parameters(k).Value   = mIn.Parameters(k).Value;
    mOut.Parameters(k).Minimum = mIn.Parameters(k).Minimum;
    mOut.Parameters(k).Maximum = mIn.Parameters(k).Maximum;
    mOut.Parameters(k).Fixed   = mIn.Parameters(k).Fixed;
end

mOut.Name = mIn.Name;

end

function m = fixPassiveInitialStates(m)
%FIXPASSIVEINITIALSTATES Configure initial states during parameter estimation.
%
% Keep theta(0) fixed to the measured first sample, but estimate theta_dot(0).
% A small nonzero theta_dot(0) mainly corrects phase/lead-lag errors caused by
% cutting the data slightly before or after a true peak.

fixedTheta = true(size(m.InitialStates(1).Value));
fixedThetaDot = false(size(m.InitialStates(2).Value));

m = setinit(m, 'Fixed', {fixedTheta; fixedThetaDot});

end

function value = getParValue(m, parName)
%GETPARVALUE Get parameter value by name.

for k = 1:length(m.Parameters)
    if strcmp(m.Parameters(k).Name, parName)
        value = m.Parameters(k).Value;
        return;
    end
end

error("Parameter '%s' not found.", parName);

end

function m = setParFixed(m, parName, fixedValue)
%SETPARFIXED Set Fixed property of parameter by name.

found = false;

for k = 1:length(m.Parameters)
    if strcmp(m.Parameters(k).Name, parName)
        m.Parameters(k).Fixed = fixedValue;
        found = true;
        break;
    end
end

if ~found
    error("Parameter '%s' not found.", parName);
end

end

function m = randomizeFreeParameters(mBase, spread)
%RANDOMIZEFREEPARAMETERS Randomize all non-fixed parameters within bounds.
%
% spread controls the log-normal perturbation strength. A value around
% 0.4--0.8 usually gives enough variation without making most starts fail.

m = mBase;

for k = 1:length(m.Parameters)

    if isParameterFixed(m.Parameters(k))
        continue;
    end

    p = m.Parameters(k).Value;
    pMin = m.Parameters(k).Minimum;
    pMax = m.Parameters(k).Maximum;

    if isempty(pMin) || isinf(pMin)
        pMin = 0;
    end

    if isempty(pMax) || isinf(pMax)
        pMax = 10 * max(abs(p), 1);
    end

    if p > 0
        pNew = p * exp(spread * randn);
    else
        pNew = pMin + rand * (pMax - pMin);
    end

    if pNew <= pMin || isnan(pNew) || isinf(pNew)
        pNew = pMin + rand * (pMax - pMin);
    end

    pNew = max(pMin, min(pNew, pMax));

    m.Parameters(k).Value = pNew;
end

end

function m = randomizeStribeckStart(mBase)
%RANDOMIZESTRIBECKSTART Randomize free Stribeck parameters within bounds.

m = randomizeFreeParameters(mBase, 0.80);

end

function isFixed = isParameterFixed(parameter)
%ISPARAMETERFIXED Robustly read idnlgrey parameter Fixed flag.

isFixed = parameter.Fixed;

if iscell(isFixed)
    isFixed = isFixed{1};
end

isFixed = all(isFixed(:));

end


function tbl = parameterUncertaintyTable(m, modelLabel)
%PARAMETERUNCERTAINTYTABLE Build a table with variance, standard deviation and 95% CI.
%
% The covariance is only meaningful for parameters that were free during
% the estimation that produced this model. Fixed parameters are reported
% with NaN variance/std because they were not statistically estimated in
% that stage.

nPar = length(m.Parameters);
model = repmat({modelLabel}, nPar, 1);
parameter = cell(nPar, 1);
value = zeros(nPar, 1);
fixed = false(nPar, 1);
variance = NaN(nPar, 1);
stdDev = NaN(nPar, 1);
relativeStdPercent = NaN(nPar, 1);
ci95Lower = NaN(nPar, 1);
ci95Upper = NaN(nPar, 1);

for k = 1:nPar
    parameter{k} = m.Parameters(k).Name;
    value(k) = double(m.Parameters(k).Value);
    fixed(k) = isParameterFixed(m.Parameters(k));
end

freeIdx = find(~fixed);

try
    covFree = getcov(m, 'value', 'free');

    if ~isempty(covFree) && isnumeric(covFree)
        if size(covFree, 1) == length(freeIdx)
            freeVariance = diag(covFree);
            variance(freeIdx) = freeVariance;
        elseif size(covFree, 1) == nPar
            variance = diag(covFree);
        end
    end
catch ME
    warning("Could not extract free-parameter covariance for model '%s': %s", modelLabel, ME.message);

    try
        covAll = getcov(m, 'value');
        if ~isempty(covAll) && isnumeric(covAll) && size(covAll, 1) == nPar
            variance = diag(covAll);
            variance(fixed) = NaN;
        end
    catch ME2
        warning("Could not extract full parameter covariance for model '%s': %s", modelLabel, ME2.message);
    end
end

% Numerical covariance calculations may produce tiny negative values on the
% diagonal due to round-off. These are clipped to zero before sqrt.
validVariance = isfinite(variance);
variance(validVariance) = max(variance(validVariance), 0);
stdDev(validVariance) = sqrt(variance(validVariance));

for k = 1:nPar
    if isfinite(stdDev(k))
        ci95Lower(k) = value(k) - 1.96 * stdDev(k);
        ci95Upper(k) = value(k) + 1.96 * stdDev(k);

        if abs(value(k)) > eps
            relativeStdPercent(k) = 100 * stdDev(k) / abs(value(k));
        end
    end
end

tbl = table(model, parameter, value, fixed, variance, stdDev, ...
    relativeStdPercent, ci95Lower, ci95Upper, ...
    'VariableNames', {'Model', 'Parameter', 'Value', 'Fixed', 'Variance', ...
    'StdDev', 'RelativeStdPercent', 'CI95Lower', 'CI95Upper'});

end

function summary = parameterIdentifiabilitySummary(m, modelLabel)
%PARAMETERIDENTIFIABILITYSUMMARY Compact identifiability diagnostics.

uncertainty = parameterUncertaintyTable(m, modelLabel);
freeRows = ~uncertainty.Fixed;

numParameters = height(uncertainty);
numFreeParameters = sum(freeRows);

if any(freeRows)
    maxRelativeStdPercent = max(uncertainty.RelativeStdPercent(freeRows), [], 'omitnan');
else
    maxRelativeStdPercent = NaN;
end

[corrFree, ~] = freeParameterCorrelationMatrix(m);

if isempty(corrFree) || size(corrFree, 1) < 2
    maxAbsCorrelation = NaN;
else
    mask = ~eye(size(corrFree));
    maxAbsCorrelation = max(abs(corrFree(mask)), [], 'omitnan');
end

try
    covFree = getcov(m, 'value', 'free');
    if isempty(covFree) || ~isnumeric(covFree) || size(covFree, 1) < 2
        covarianceConditionNumber = NaN;
    else
        covarianceConditionNumber = cond(covFree);
    end
catch
    covarianceConditionNumber = NaN;
end

summary = table({modelLabel}, numParameters, numFreeParameters, ...
    maxRelativeStdPercent, maxAbsCorrelation, covarianceConditionNumber, ...
    'VariableNames', {'Model', 'NumParameters', 'NumFreeParameters', ...
    'MaxRelativeStdPercent', 'MaxAbsFreeParameterCorrelation', ...
    'FreeCovarianceConditionNumber'});

end

function corrTbl = parameterCorrelationTable(m)
%PARAMETERCORRELATIONTABLE Return the free-parameter correlation matrix as a table.

[corrFree, freeNames] = freeParameterCorrelationMatrix(m);

if isempty(corrFree)
    corrTbl = table();
    return;
end

rowNames = matlab.lang.makeUniqueStrings(freeNames);
varNames = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(freeNames));

corrTbl = array2table(corrFree, 'VariableNames', varNames, 'RowNames', rowNames);

end

function [corrFree, freeNames] = freeParameterCorrelationMatrix(m)
%FREEPARAMETERCORRELATIONMATRIX Compute correlation matrix of free parameters.

nPar = length(m.Parameters);
fixed = false(nPar, 1);
parameterNames = cell(nPar, 1);

for k = 1:nPar
    fixed(k) = isParameterFixed(m.Parameters(k));
    parameterNames{k} = m.Parameters(k).Name;
end

freeIdx = find(~fixed);
freeNames = parameterNames(freeIdx);

if isempty(freeIdx)
    corrFree = [];
    freeNames = {};
    return;
end

try
    covFree = getcov(m, 'value', 'free');
catch
    corrFree = [];
    freeNames = {};
    return;
end

if isempty(covFree) || ~isnumeric(covFree)
    corrFree = [];
    freeNames = {};
    return;
end

if size(covFree, 1) ~= length(freeIdx)
    if size(covFree, 1) == nPar
        covFree = covFree(freeIdx, freeIdx);
    else
        corrFree = [];
        freeNames = {};
        return;
    end
end

stdFree = sqrt(max(diag(covFree), 0));
denominator = stdFree * stdFree.';

corrFree = covFree ./ denominator;
corrFree(denominator == 0) = NaN;
corrFree(1:size(corrFree,1)+1:end) = 1;

end

function dispPassiveParameters(m)
%DISPPASSIVEPARAMETERS Display parameter values.

fprintf("\nEstimated parameters for model: %s\n", m.Name);

for k = 1:length(m.Parameters)
    fprintf("  %-12s = %.8g", m.Parameters(k).Name, m.Parameters(k).Value);

    if isParameterFixed(m.Parameters(k))
        fprintf("  fixed\n");
    else
        fprintf("  free\n");
    end
end

end

function [yHatData, fit] = comparePassiveModel(data, model)
%COMPAREPASSIVEMODEL Compare using estimated initial conditions when available.
%
% This is useful for validation and plotting because the manually selected
% start index can be slightly before/after a peak. If compareOptions with
% InitialCondition='estimate' is unavailable in the MATLAB version, fall back
% to the standard compare call.

try
    optCompare = compareOptions;
    optCompare.InitialCondition = 'estimate';
    [yHatData, fit] = compare(data, model, optCompare);
catch
    try
        optCompare = compareOptions('InitialCondition', 'estimate');
        [yHatData, fit] = compare(data, model, optCompare);
    catch
        [yHatData, fit] = compare(data, model);
    end
end

end

function m = meanFitValue(fit)
%MEANFITVALUE Convert compare(...) fit output to scalar mean.

values = collectFitValues(fit);
values = values(isfinite(values));

if isempty(values)
    m = NaN;
else
    m = mean(values);
end

end

function values = collectFitValues(fit)
%COLLECTFITVALUES Recursively collect numeric fit values.

if isnumeric(fit)
    values = fit(:);
elseif iscell(fit)
    values = [];

    for k = 1:numel(fit)
        values = [values; collectFitValues(fit{k})]; %#ok<AGROW>
    end
else
    error("Unsupported fit output type: %s", class(fit));
end

end

function plotModelAgainstAllMeasurements(model, modelName, zCell, expData, angles, segment, showTailStart)
%PLOTMODELAGAINSTALLMEASUREMENTS Plot measured and simulated data in 4x2 layout.

figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(angles)

    nexttile;

    mSingle = setPassiveInitialStates(model, expData, i, segment);

    try
        [yHatData, fit] = comparePassiveModel(zCell{i}, mSingle);
        yModel = iddataOutputVector(yHatData);
        fitValue = meanFitValue(fit);
    catch ME
        warning("Could not simulate model '%s' for experiment %d: %s", modelName, i, ME.message);
        yModel = NaN;
        fitValue = NaN;
    end

    switch segment
        case 'id'
            t = expData(i).t_id(:);
            y = expData(i).y_id(:);
        case 'tail'
            t = expData(i).t_tail(:);
            y = expData(i).y_tail(:);
        otherwise
            error("Unknown segment '%s'.", segment);
    end

    plot(t, y, 'DisplayName', 'Measured');
    hold on;

    if numel(yModel) > 1
        n = min(length(t), length(yModel));
        plot(t(1:n), yModel(1:n), '--', 'LineWidth', 1.1, 'DisplayName', 'Model');
    end

    if showTailStart && strcmp(segment, 'id')
        xline(expData(i).tail_start_time, ':', 'Tail start', 'HandleVisibility', 'off');
    end

    title(sprintf('%d deg, fit %.1f%%', angles(i), fitValue));
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('Location', 'best');
end

sgtitle(modelName);

end

function y = iddataOutputVector(dataObject)
%IDDATAOUTPUTVECTOR Extract output data from iddata returned by compare.

try
    y = get(dataObject, 'OutputData');
catch
    y = dataObject.OutputData;
end

if iscell(y)
    y = y{1};
end

y = double(y(:));

end

function writePassiveModelFiles(folder, overwriteFiles)
%WRITEPASSIVEMODELFILES Write idnlgrey model files used by this script.

writeTextFile(fullfile(folder, 'passive_rod_viscous_model.m'), passiveViscousModelText(), overwriteFiles);
writeTextFile(fullfile(folder, 'passive_rod_coulomb_model.m'), passiveCoulombModelText(), overwriteFiles);
writeTextFile(fullfile(folder, 'passive_rod_stribeck_report_model.m'), passiveStribeckReportModelText(), overwriteFiles);

end

function writeTextFile(fileName, lines, overwriteFile)
%WRITETEXTFILE Write a cell array of text lines to file.

if exist(fileName, 'file') && ~overwriteFile
    return;
end

fid = fopen(fileName, 'w');

if fid < 0
    error("Could not write file: %s", fileName);
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

for k = 1:numel(lines)
    fprintf(fid, '%s\n', lines{k});
end

end

function lines = passiveViscousModelText()
%PASSIVEVISCOUSMODELTEXT Text for passive_rod_viscous_model.m.

lines = {
'function [dx, y] = passive_rod_viscous_model(t, x, u, p_b2, p_g2, varargin)'
'%PASSIVE_ROD_VISCOUS_MODEL Passive rod model with viscous damping.'
'% theta_ddot = -p_b2 theta_dot - p_g2 sin(theta)'
''
'psi = x(1);'
'theta_dot = x(2);'
''
'dx = zeros(2, 1);'
'dx(1) = theta_dot;'
'dx(2) = -p_b2 * theta_dot - p_g2 * sin(psi);'
''
'y = psi;'
''
'end'
};

end

function lines = passiveCoulombModelText()
%PASSIVECOULOMBMODELTEXT Text for passive_rod_coulomb_model.m.

lines = {
'function [dx, y] = passive_rod_coulomb_model(t, x, u, p_b2, p_g2, p_c2, eps_v, varargin)'
'%PASSIVE_ROD_COULOMB_MODEL Passive rod model with Coulomb-viscous friction.'
'% theta_ddot = -p_b2 theta_dot - p_g2 sin(theta) - p_c2 tanh(theta_dot/eps_v)'
''
'psi = x(1);'
'theta_dot = x(2);'
''
'eps_v = max(eps_v, 1e-5);'
''
'dx = zeros(2, 1);'
'dx(1) = theta_dot;'
'dx(2) = -p_b2 * theta_dot ...'
'        -p_g2 * sin(psi) ...'
'        -p_c2 * tanh(theta_dot / eps_v);'
''
'y = psi;'
''
'end'
};

end

function lines = passiveStribeckReportModelText()
%PASSIVESTRIBECKREPORTMODELTEXT Text for passive_rod_stribeck_report_model.m.

lines = {
'function [dx, y] = passive_rod_stribeck_report_model(t, x, u, p_b2, p_g2, p_c2, p_sdelta, v_s, eps_v, varargin)'
'%PASSIVE_ROD_STRIBECK_REPORT_MODEL Report-derived passive rod Stribeck model.'
'% p_sdelta = F_S - F_C, so the static level is p_c2 + p_sdelta.'
'% F_Str = p_b2*theta_dot + (p_c2 + p_sdelta*exp(-(theta_dot/v_s)^2))*tanh(theta_dot/eps_v)'
'% theta_ddot = -F_Str - p_g2*sin(theta)'
''
'psi = x(1);'
'theta_dot = x(2);'
''
'v_s = max(v_s, 1e-5);'
'eps_v = max(eps_v, 1e-5);'
''
'friction_level = p_c2 + p_sdelta * exp(-(theta_dot / v_s)^2);'
'friction = p_b2 * theta_dot + friction_level * tanh(theta_dot / eps_v);'
''
'dx = zeros(2, 1);'
'dx(1) = theta_dot;'
'dx(2) = -friction - p_g2 * sin(psi);'
''
'y = psi;'
''
'end'
};

end
