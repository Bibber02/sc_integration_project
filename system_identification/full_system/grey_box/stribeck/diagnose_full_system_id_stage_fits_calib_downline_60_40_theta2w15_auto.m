clear;
clc;
close all;

%% ================================================================
% Diagnostics for full-system grey-box identification results
%
% This script is separate from the estimation script. It loads the saved
% stage result .mat files, recreates correctly sized one-experiment models,
% compares them to every training and validation experiment, plots tiled
% measured-vs-model figures, and writes a diagnostic table with:
%
%   - MATLAB fit percentage
%   - VAF, RMSE, and NRMSE as used in the ID lecture
%   - RMSE in rad and deg
%   - mean residual / bias in rad and deg
%   - RMSE after removing mean residual
%   - amplitude ratio model/measured
%   - phase lag from cross-correlation
%   - residual/input correlation
%
% The comparison uses InitialCondition = 'model'. This is stricter than
% estimating the initial condition during plotting.
% ================================================================

%% Settings

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end

modelsFolder = scriptFolder;
addpath(modelsFolder);

modelFile = 'greybox_id_full_stribeck_model_calib_downline_no_p0';

% Result files produced by the AUTO-search full-system script.
stage2MatFile = fullfile(scriptFolder, ...
    'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage2_active_offsets.mat');
stage3MatFile = fullfile(scriptFolder, ...
    'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_sensor_calibration.mat');
stage4MatFile = fullfile(scriptFolder, ...
    'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage4_final_refinement.mat');

% If you want to analyse the older files instead, change the filenames above.

diagnosticsCsvFile = fullfile(scriptFolder, ...
    'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_diagnostics.csv');
summaryCsvFile = fullfile(scriptFolder, ...
    'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_diagnostics_summary.csv');

makePlots = true;
figureUnits = 'pixels';
figurePosition = [60 60 1450 850];

maxLagSeconds = 1.0;

compareOpt = compareOptions;
compareOpt.InitialCondition = 'model';

%% Load saved result files into separate structs

S2 = load(stage2MatFile);
S3 = load(stage3MatFile);
hasStage4 = isfile(stage4MatFile);
if hasStage4
    S4 = load(stage4MatFile);
end

zEst = S2.zEst;
zVal = S2.zVal;
x0Est = S2.x0Est;
x0Val = S2.x0Val;
experimentNamesEst = S2.experimentNamesEst;
experimentNamesVal = S2.experimentNamesVal;

order = S2.order;
TsModel = S2.TsModel;
parameterNames = S2.parameterNames;
minimumValues = S2.minimumValues;
maximumValues = S2.maximumValues;

%% Define stages to evaluate

stageNames = {'Stage 1', 'Stage 2', 'Stage 3'};
stageModels = {S2.modelStage1, S2.modelStage2, S3.modelStage3};

if hasStage4
    stageNames{end+1} = 'Stage 4'; %#ok<SAGROW>
    stageModels{end+1} = S4.modelStage4; %#ok<SAGROW>
end

%% Evaluate and plot all stages

allDiagnostics = table();

for kStage = 1:numel(stageNames)
    stageName = stageNames{kStage};
    stageModel = stageModels{kStage};

    fprintf('\n======================================================\n');
    fprintf('%s diagnostics\n', stageName);
    fprintf('======================================================\n');

    [fitTrain, diagTrain] = evaluateStageOnDataset(zEst, stageModel, x0Est, experimentNamesEst, ...
        modelFile, order, TsModel, parameterNames, minimumValues, maximumValues, compareOpt, ...
        stageName, 'training', makePlots, figureUnits, figurePosition, maxLagSeconds);

    [fitVal, diagVal] = evaluateStageOnDataset(zVal, stageModel, x0Val, experimentNamesVal, ...
        modelFile, order, TsModel, parameterNames, minimumValues, maximumValues, compareOpt, ...
        stageName, 'validation', makePlots, figureUnits, figurePosition, maxLagSeconds);

    allDiagnostics = [allDiagnostics; diagTrain; diagVal]; %#ok<AGROW>

    fprintf('%-20s training mean fit   = %8.3f %%\n', stageName, mean(flattenFitValues(fitTrain)));
    fprintf('%-20s validation mean fit = %8.3f %%\n', stageName, mean(flattenFitValues(fitVal)));
    printCalibrationSummary(stageName, stageModel);
end

%% Save detailed diagnostics

writetable(allDiagnostics, diagnosticsCsvFile);
fprintf('\nSaved detailed diagnostics:\n  %s\n', diagnosticsCsvFile);

%% Build and save compact summary table

summaryRows = {};
summaryStage = {};
summarySet = {};
summaryOutput = {};
summaryMeanFit = [];
summaryMeanVAF = [];
summaryMeanNRMSE = [];
summaryMeanNRMSEFit = [];
summaryMeanRmseDeg = [];
summaryMeanBiasDeg = [];
summaryMeanRmseNoBiasDeg = [];
summaryMeanAmpRatio = [];
summaryMeanAbsPhaseLag = [];
summaryMeanAbsResidInputCorr = [];

uniqueStages = unique(allDiagnostics.Stage, 'stable');
uniqueSets = unique(allDiagnostics.DataSet, 'stable');
uniqueOutputs = unique(allDiagnostics.Output, 'stable');

for iStage = 1:numel(uniqueStages)
    for iSet = 1:numel(uniqueSets)
        for iOut = 1:numel(uniqueOutputs)
            mask = strcmp(allDiagnostics.Stage, uniqueStages{iStage}) & ...
                   strcmp(allDiagnostics.DataSet, uniqueSets{iSet}) & ...
                   strcmp(allDiagnostics.Output, uniqueOutputs{iOut});
            if ~any(mask)
                continue;
            end

            summaryStage{end+1,1} = uniqueStages{iStage}; %#ok<SAGROW>
            summarySet{end+1,1} = uniqueSets{iSet}; %#ok<SAGROW>
            summaryOutput{end+1,1} = uniqueOutputs{iOut}; %#ok<SAGROW>
            summaryMeanFit(end+1,1) = mean(allDiagnostics.FitPercent(mask), 'omitnan'); %#ok<SAGROW>
            summaryMeanVAF(end+1,1) = mean(allDiagnostics.VAF_Percent(mask), 'omitnan'); %#ok<SAGROW>
            summaryMeanNRMSE(end+1,1) = mean(allDiagnostics.NRMSE(mask), 'omitnan'); %#ok<SAGROW>
            summaryMeanNRMSEFit(end+1,1) = mean(allDiagnostics.NRMSE_FitPercent(mask), 'omitnan'); %#ok<SAGROW>
            summaryMeanRmseDeg(end+1,1) = mean(allDiagnostics.RMSE_deg(mask), 'omitnan'); %#ok<SAGROW>
            summaryMeanBiasDeg(end+1,1) = mean(allDiagnostics.Bias_deg(mask), 'omitnan'); %#ok<SAGROW>
            summaryMeanRmseNoBiasDeg(end+1,1) = mean(allDiagnostics.RMSE_NoBias_deg(mask), 'omitnan'); %#ok<SAGROW>
            summaryMeanAmpRatio(end+1,1) = mean(allDiagnostics.AmplitudeRatio(mask), 'omitnan'); %#ok<SAGROW>
            summaryMeanAbsPhaseLag(end+1,1) = mean(abs(allDiagnostics.PhaseLag_s(mask)), 'omitnan'); %#ok<SAGROW>
            summaryMeanAbsResidInputCorr(end+1,1) = mean(abs(allDiagnostics.ResidInputCorrMaxAbs(mask)), 'omitnan'); %#ok<SAGROW>
        end
    end
end

summaryTable = table(summaryStage, summarySet, summaryOutput, summaryMeanFit, ...
    summaryMeanVAF, summaryMeanNRMSE, summaryMeanNRMSEFit, ...
    summaryMeanRmseDeg, summaryMeanBiasDeg, summaryMeanRmseNoBiasDeg, ...
    summaryMeanAmpRatio, summaryMeanAbsPhaseLag, summaryMeanAbsResidInputCorr, ...
    'VariableNames', {'Stage', 'DataSet', 'Output', 'MeanFitPercent', ...
    'MeanVAF_Percent', 'MeanNRMSE', 'MeanNRMSE_FitPercent', ...
    'MeanRMSE_deg', 'MeanBias_deg', 'MeanRMSE_NoBias_deg', ...
    'MeanAmplitudeRatio', 'MeanAbsPhaseLag_s', 'MeanAbsResidualInputCorrelation'});

fprintf('\n======================================================\n');
fprintf('Diagnostic summary\n');
fprintf('======================================================\n\n');
disp(summaryTable);

writetable(summaryTable, summaryCsvFile);
fprintf('\nSaved diagnostic summary:\n  %s\n', summaryCsvFile);

%% ================================================================
% Local functions
% ================================================================

function [fitAll, diagnostics] = evaluateStageOnDataset(zData, modelStage, x0, experimentNames, ...
    modelFile, order, TsModel, parameterNames, minimumValues, maximumValues, compareOpt, ...
    stageName, dataSetName, makePlots, figureUnits, figurePosition, maxLagSeconds)

    nExp = size(x0, 2);
    fitAll = cell(nExp, 1);

    if makePlots
        figName = sprintf('%s %s fit', stageName, dataSetName);
        figure('Name', figName, 'Units', figureUnits, 'Position', figurePosition);
        nCols = ceil(sqrt(nExp));
        nRows = ceil(nExp / nCols);
        tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    diagnostics = table();

    for kExp = 1:nExp
        zExp = getexp(zData, kExp);

        modelExp = makeModelForOneExperiment(modelStage, modelFile, order, x0(:, kExp), TsModel, ...
            parameterNames, minimumValues, maximumValues);

        [zModelExp, fitExp] = compare(zExp, modelExp, compareOpt);
        fitAll{kExp} = fitExp;

        yMeasured = getOutputArray(zExp);
        yModel = getOutputArray(zModelExp);
        u = getInputArray(zExp);

        nSamples = min([size(yMeasured, 1), size(yModel, 1), length(u)]);
        yMeasured = yMeasured(1:nSamples, :);
        yModel = yModel(1:nSamples, :);
        u = u(1:nSamples, :);

        t = getTimeVector(zExp, nSamples);
        Ts = getTsScalar(zExp);

        if makePlots
            nexttile;
            plot(t, yMeasured(:, 1), 'LineWidth', 1.0);
            hold on;
            plot(t, yModel(:, 1), '--', 'LineWidth', 1.0);
            plot(t, yMeasured(:, 2), 'LineWidth', 1.0);
            plot(t, yModel(:, 2), '--', 'LineWidth', 1.0);
            grid on;

            fitFlat = flattenFitValues(fitExp);
            if isempty(fitFlat)
                fitText = '';
            else
                fitText = sprintf(' | mean %.1f%%', mean(fitFlat));
            end

            title([experimentNames{kExp}, fitText], 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Angle [rad]');

            if kExp == 1
                legend({'\theta_1 measured', '\theta_1 model', '\theta_2 measured', '\theta_2 model'}, ...
                    'Location', 'best');
            end
        end

        fitValues = flattenFitValues(fitExp);
        if numel(fitValues) < 2
            fitValues = [NaN; NaN];
        end

        for outputIndex = 1:2
            yM = yMeasured(:, outputIndex);
            yH = yModel(:, outputIndex);
            e = yM - yH;

            rmseRad = sqrt(mean(e.^2, 'omitnan'));
            biasRad = mean(e, 'omitnan');
            eNoBias = e - biasRad;
            rmseNoBiasRad = sqrt(mean(eNoBias.^2, 'omitnan'));

            % Lecture 6 validation metrics:
            %   VAF   = max(0, (1 - ||y-yhat||_2^2 / ||y||_2^2) * 100%)
            %   RMSE  = sqrt(mean((y-yhat).^2))
            %   NRMSE = ||y-yhat||_2 / ||y-mean(y)||_2
            eFinite = e(isfinite(e));
            yMFinite = yM(isfinite(e) & isfinite(yM));
            if isempty(eFinite) || isempty(yMFinite)
                vafPercent = NaN;
                nrmseLecture = NaN;
                nrmseFitPercent = NaN;
            else
                errorNormSquared = sum(eFinite.^2);
                yNormSquared = sum(yMFinite.^2);
                yCenteredNorm = norm(yMFinite - mean(yMFinite));
                if yNormSquared > eps
                    vafPercent = max(0, (1 - errorNormSquared / yNormSquared) * 100);
                else
                    vafPercent = NaN;
                end
                if yCenteredNorm > eps
                    nrmseLecture = norm(eFinite) / yCenteredNorm;
                    nrmseFitPercent = 100 * (1 - nrmseLecture);
                else
                    nrmseLecture = NaN;
                    nrmseFitPercent = NaN;
                end
            end

            ampMeasured = 0.5 * (max(yM) - min(yM));
            ampModel = 0.5 * (max(yH) - min(yH));
            if ampMeasured > eps
                ampRatio = ampModel / ampMeasured;
            else
                ampRatio = NaN;
            end

            [phaseLagSeconds, peakCorr] = estimatePhaseLag(yM, yH, Ts, maxLagSeconds);
            [residInputCorr0, residInputCorrMaxAbs] = residualInputCorrelation(e, u(:, 1), Ts, maxLagSeconds);

            thisRow = table({stageName}, {dataSetName}, {experimentNames{kExp}}, ...
                {signalTypeFromName(experimentNames{kExp})}, kExp, {sprintf('theta_%d', outputIndex)}, ...
                fitValues(outputIndex), vafPercent, nrmseLecture, nrmseFitPercent, ...
                rmseRad, rad2deg(rmseRad), biasRad, rad2deg(biasRad), ...
                rmseNoBiasRad, rad2deg(rmseNoBiasRad), ampMeasured, ampModel, ampRatio, ...
                phaseLagSeconds, peakCorr, residInputCorr0, residInputCorrMaxAbs, ...
                'VariableNames', {'Stage', 'DataSet', 'Experiment', 'SignalType', ...
                'ExperimentIndex', 'Output', 'FitPercent', 'VAF_Percent', 'NRMSE', 'NRMSE_FitPercent', ...
                'RMSE_rad', 'RMSE_deg', 'Bias_rad', 'Bias_deg', ...
                'RMSE_NoBias_rad', 'RMSE_NoBias_deg', ...
                'AmplitudeMeasured_rad', 'AmplitudeModel_rad', 'AmplitudeRatio', ...
                'PhaseLag_s', 'PhasePeakCorrelation', 'ResidInputCorr0', 'ResidInputCorrMaxAbs'});

            diagnostics = [diagnostics; thisRow]; %#ok<AGROW>
        end
    end

    if makePlots
        sgtitle(sprintf('%s %s fit', stageName, dataSetName));
    end
end

function modelOut = makeModelForOneExperiment(modelIn, modelFile, order, x0, TsModel, ...
    parameterNames, minimumValues, maximumValues)

    parameters = cell(length(modelIn.Parameters), 1);
    for kPar = 1:length(modelIn.Parameters)
        parameters{kPar} = modelIn.Parameters(kPar).Value;
    end

    initialStates = {x0(1); x0(2); x0(3); x0(4)};
    modelOut = idnlgrey(modelFile, order, parameters, initialStates, TsModel);

    for kPar = 1:length(parameterNames)
        modelOut.Parameters(kPar).Name = parameterNames{kPar};
        modelOut.Parameters(kPar).Minimum = minimumValues(kPar);
        modelOut.Parameters(kPar).Maximum = maximumValues(kPar);
        modelOut.Parameters(kPar).Fixed = true;
    end

    modelOut.InputName = modelIn.InputName;
    modelOut.InputUnit = modelIn.InputUnit;
    modelOut.OutputName = modelIn.OutputName;
    modelOut.OutputUnit = modelIn.OutputUnit;
    modelOut.TimeUnit = modelIn.TimeUnit;

    for kState = 1:4
        modelOut.InitialStates(kState).Name = modelIn.InitialStates(kState).Name;
        modelOut.InitialStates(kState).Unit = modelIn.InitialStates(kState).Unit;
    end

    modelOut = setinit(modelOut, 'Fixed', {true; true; true; true});
end

function y = getOutputArray(z)
    if iscell(z)
        z = z{1};
    end

    if isa(z, 'iddata')
        y = get(z, 'OutputData');
        if iscell(y), y = y{1}; end
    elseif isnumeric(z)
        y = z;
    else
        try
            y = get(z, 'OutputData');
            if iscell(y), y = y{1}; end
        catch
            y = z.OutputData;
            if iscell(y), y = y{1}; end
        end
    end

    y = squeeze(y);
    if isrow(y) && numel(y) > 2
        y = y(:);
    end
end

function u = getInputArray(z)
    if isa(z, 'iddata')
        u = get(z, 'InputData');
        if iscell(u), u = u{1}; end
    else
        u = get(z, 'InputData');
        if iscell(u), u = u{1}; end
    end
    u = squeeze(u);
    if isrow(u)
        u = u(:);
    end
end

function t = getTimeVector(z, nSamples)
    t = [];
    try
        samplingInstants = get(z, 'SamplingInstants');
        if iscell(samplingInstants), samplingInstants = samplingInstants{1}; end
        if numel(samplingInstants) >= nSamples
            t = samplingInstants(1:nSamples);
            t = t(:);
        end
    catch
    end
    if isempty(t)
        Ts = getTsScalar(z);
        t = (0:nSamples-1)' * Ts;
    end
end

function Ts = getTsScalar(z)
    Ts = get(z, 'Ts');
    if iscell(Ts), Ts = Ts{1}; end
end

function [lagSeconds, peakCorr] = estimatePhaseLag(yMeasured, yModel, Ts, maxLagSeconds)
    y1 = yMeasured(:) - mean(yMeasured(:), 'omitnan');
    y2 = yModel(:) - mean(yModel(:), 'omitnan');

    if std(y1, 'omitnan') < eps || std(y2, 'omitnan') < eps
        lagSeconds = NaN;
        peakCorr = NaN;
        return;
    end

    maxLagSamples = max(1, round(maxLagSeconds / Ts));
    [c, lags] = xcorr(y1, y2, maxLagSamples, 'coeff');
    [peakCorr, idx] = max(c);
    lagSeconds = lags(idx) * Ts;
end

function [corr0, corrMaxAbs] = residualInputCorrelation(e, u, Ts, maxLagSeconds)
    e = e(:) - mean(e(:), 'omitnan');
    u = u(:) - mean(u(:), 'omitnan');

    if std(e, 'omitnan') < eps || std(u, 'omitnan') < eps
        corr0 = NaN;
        corrMaxAbs = NaN;
        return;
    end

    C = corrcoef(e, u);
    corr0 = C(1, 2);

    maxLagSamples = max(1, round(maxLagSeconds / Ts));
    c = xcorr(e, u, maxLagSamples, 'coeff');
    corrMaxAbs = max(abs(c));
end

function values = flattenFitValues(fitValue)
    values = [];
    if iscell(fitValue)
        for k = 1:numel(fitValue)
            values = [values; flattenFitValues(fitValue{k})]; %#ok<AGROW>
        end
    else
        values = fitValue(:);
    end
    values = values(isfinite(values));
end

function signalType = signalTypeFromName(experimentName)
    if contains(lower(experimentName), 'prbs')
        signalType = 'PRBS';
    elseif contains(lower(experimentName), 'chirp')
        signalType = 'chirp';
    else
        signalType = 'unknown';
    end
end

function printCalibrationSummary(stageName, modelIn)
    thetaScale = getParameterValue(modelIn, 'theta_scale');
    theta1Offset = getParameterValue(modelIn, 'theta1_offset');
    thetaAbsDownRaw = getParameterValue(modelIn, 'theta_abs_down_raw');
    theta2Offset = pi - thetaScale*thetaAbsDownRaw - theta1Offset;

    fprintf('%s theta_scale  = %.9g\n', stageName, thetaScale);
    fprintf('%s theta1_offset = %.9g rad = %.6g deg\n', stageName, theta1Offset, rad2deg(theta1Offset));
    fprintf('%s theta_abs_down_raw = %.9g rad = %.6g deg\n', stageName, thetaAbsDownRaw, rad2deg(thetaAbsDownRaw));
    fprintf('%s theta2_offset = %.9g rad = %.6g deg   computed from down-line constraint\n', stageName, theta2Offset, rad2deg(theta2Offset));
    fprintf('%s down-line check at raw theta1+theta2 = theta_abs_down_raw: %.9g rad = %.6g deg\n', ...
        stageName, thetaScale*thetaAbsDownRaw + theta1Offset + theta2Offset, rad2deg(thetaScale*thetaAbsDownRaw + theta1Offset + theta2Offset));
end

function value = getParameterValue(modelIn, parameterName)
    value = NaN;
    for kPar = 1:length(modelIn.Parameters)
        if strcmp(char(modelIn.Parameters(kPar).Name), parameterName)
            value = modelIn.Parameters(kPar).Value;
            return;
        end
    end
end
