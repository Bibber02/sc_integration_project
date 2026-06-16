clear;
clc;
close all;

%% ================================================================
% Diagnose initial theta_2 oscillation in full-system grey-box models
%
% Purpose:
%   This script checks whether the early theta_2 oscillation is caused by
%   a nonzero model acceleration at the measured initial state, even though
%   the initial velocities are set to zero.
%
% It loads the saved Stage 2/3/4 workspaces from the calibration ID script,
% evaluates dx/dt at t = 0 for every experiment, and compares the measured
% and simulated theta_2 motion over the first few seconds.
%
% Run this script from:
%   system_identification/full_system/grey_box/stribeck
%% ================================================================

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end
addpath(scriptFolder);

modelFile = 'greybox_id_full_stribeck_model_calib_downline_no_p0';
order = [2 1 4];
TsModel = 0;

stageFiles = {
    'Stage 2', fullfile(scriptFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage2_active_offsets.mat'), 'modelStage2';
    'Stage 3', fullfile(scriptFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage3_sensor_calibration.mat'), 'modelStage3';
    'Stage 4', fullfile(scriptFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage4_final_refinement.mat'), 'modelStage4'};

zoomSeconds = 6.0;
restSecondsForMetrics = 4.0;
makeZoomPlots = true;
stagesToPlot = {'Stage 3', 'Stage 4'};
dataSetsToPlot = {'validation'};

compareOpt = compareOptions;
compareOpt.InitialCondition = 'model';

allRows = table();

fprintf('\n======================================================\n');
fprintf('Initial theta_2 oscillation diagnostic\n');
fprintf('======================================================\n\n');

whichModel = which(modelFile);
if isempty(whichModel)
    error('Cannot find model file %s on the MATLAB path.', modelFile);
end
fprintf('Using model file:\n  %s\n\n', whichModel);

for kStage = 1:size(stageFiles, 1)
    stageName = stageFiles{kStage, 1};
    stageMatFile = stageFiles{kStage, 2};
    modelVarName = stageFiles{kStage, 3};

    if ~exist(stageMatFile, 'file')
        warning('Skipping %s because result file was not found:\n  %s', stageName, stageMatFile);
        continue;
    end

    S = load(stageMatFile);

    if ~isfield(S, modelVarName)
        warning('Skipping %s because variable %s was not found in:\n  %s', stageName, modelVarName, stageMatFile);
        continue;
    end

    model = S.(modelVarName);
    parameterValues = getParameterValues(model);

    fprintf('======================================================\n');
    fprintf('%s\n', stageName);
    fprintf('======================================================\n');
    thetaScaleStage = getParameterValueByName(model, 'theta_scale');
    theta1OffsetStage = getParameterValueByName(model, 'theta1_offset');
    thetaAbsDownRawStage = getParameterValueByName(model, 'theta_abs_down_raw');
    theta2OffsetStage = pi - thetaScaleStage*thetaAbsDownRawStage - theta1OffsetStage;

    fprintf('theta_scale  = %.9g\n', thetaScaleStage);
    fprintf('theta1_offset = %.9g rad = %.6g deg\n', theta1OffsetStage, rad2deg(theta1OffsetStage));
    fprintf('theta_abs_down_raw = %.9g rad = %.6g deg\n', thetaAbsDownRawStage, rad2deg(thetaAbsDownRawStage));
    fprintf('theta2_offset = %.9g rad = %.6g deg   computed from down-line constraint\n', theta2OffsetStage, rad2deg(theta2OffsetStage));
    fprintf('down-line check at raw theta1+theta2 = theta_abs_down_raw: %.9g rad = %.6g deg\n\n', ...
        thetaScaleStage*thetaAbsDownRawStage + theta1OffsetStage + theta2OffsetStage, ...
        rad2deg(thetaScaleStage*thetaAbsDownRawStage + theta1OffsetStage + theta2OffsetStage));

    dataSets = {
        'training',   S.zEst, S.x0Est, S.experimentNamesEst;
        'validation', S.zVal, S.x0Val, S.experimentNamesVal};

    for kSet = 1:size(dataSets, 1)
        dataSetName = dataSets{kSet, 1};
        zData = dataSets{kSet, 2};
        x0Data = dataSets{kSet, 3};
        experimentNames = dataSets{kSet, 4};
        nExp = size(x0Data, 2);

        theta2DdotValues = NaN(nExp, 1);
        theta2ModelRangeValues = NaN(nExp, 1);
        theta2MeasuredRangeValues = NaN(nExp, 1);

        if makeZoomPlots && any(strcmp(stageName, stagesToPlot)) && any(strcmp(dataSetName, dataSetsToPlot))
            figure('Name', sprintf('%s %s initial theta_2 zoom', stageName, dataSetName), ...
                'Units', 'pixels', 'Position', [80 80 1500 850]);
            nCols = ceil(sqrt(nExp));
            nRows = ceil(nExp / nCols);
            tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
        end

        for kExp = 1:nExp
            zExp = getExperiment(zData, kExp);
            x0 = x0Data(:, kExp);
            uArray = getInputArray(zExp);
            yMeasured = getOutputArray(zExp);
            t = getTimeVector(zExp, size(yMeasured, 1));

            u0 = uArray(1, :).';
            [dx0, ~] = feval(modelFile, 0, x0, u0, parameterValues{:});

            thetaScale = getParameterValueByName(model, 'theta_scale');
            theta1Offset = getParameterValueByName(model, 'theta1_offset');
            thetaAbsDownRaw = getParameterValueByName(model, 'theta_abs_down_raw');
            theta2Offset = pi - thetaScale*thetaAbsDownRaw - theta1Offset;

            theta1Phys0 = thetaScale * x0(1) + theta1Offset;
            theta2Phys0 = thetaScale * x0(2) + theta2Offset;
            thetaSumPhys0 = theta1Phys0 + theta2Phys0;

            modelOne = makeSingleExperimentModel(modelFile, order, parameterValues, x0, TsModel, model);
            [zModelExp, fitExp] = compare(zExp, modelOne, compareOpt); %#ok<ASGLU>
            yModel = getOutputArray(zModelExp);

            idxRest = t <= restSecondsForMetrics;
            idxZoom = t <= zoomSeconds;
            if ~any(idxRest)
                idxRest = 1:min(size(yMeasured, 1), round(restSecondsForMetrics / mean(diff(t))) + 1);
            end

            theta2MeasuredRest = yMeasured(idxRest, 2);
            theta2ModelRest = yModel(idxRest, 2);
            theta2MeasuredRange = max(theta2MeasuredRest) - min(theta2MeasuredRest);
            theta2ModelRange = max(theta2ModelRest) - min(theta2ModelRest);
            theta2InitialModelDrift1s = 0.5 * dx0(4) * 1.0^2;
            theta2InitialModelDrift025s = 0.5 * dx0(4) * 0.25^2;

            expName = experimentNames{kExp};
            newRow = table({stageName}, {dataSetName}, {expName}, ...
                x0(1), x0(2), u0(1), ...
                theta1Phys0, theta2Phys0, thetaSumPhys0, ...
                sin(theta1Phys0), sin(thetaSumPhys0), ...
                dx0(3), dx0(4), ...
                rad2deg(dx0(3)), rad2deg(dx0(4)), ...
                theta2InitialModelDrift025s, rad2deg(theta2InitialModelDrift025s), ...
                theta2InitialModelDrift1s, rad2deg(theta2InitialModelDrift1s), ...
                theta2MeasuredRange, rad2deg(theta2MeasuredRange), ...
                theta2ModelRange, rad2deg(theta2ModelRange), ...
                'VariableNames', {'Stage', 'DataSet', 'Experiment', ...
                'Theta1Raw0_rad', 'Theta2Raw0_rad', 'Input0', ...
                'Theta1Phys0_rad', 'Theta2Phys0_rad', 'Theta1PlusTheta2Phys0_rad', ...
                'SinTheta1Phys0', 'SinTheta1PlusTheta2Phys0', ...
                'Theta1Ddot0_rad_s2', 'Theta2Ddot0_rad_s2', ...
                'Theta1Ddot0_deg_s2', 'Theta2Ddot0_deg_s2', ...
                'Theta2ConstantAccelDrift025s_rad', 'Theta2ConstantAccelDrift025s_deg', ...
                'Theta2ConstantAccelDrift1s_rad', 'Theta2ConstantAccelDrift1s_deg', ...
                'Theta2MeasuredRangeFirstRest_rad', 'Theta2MeasuredRangeFirstRest_deg', ...
                'Theta2ModelRangeFirstRest_rad', 'Theta2ModelRangeFirstRest_deg'});

            allRows = [allRows; newRow]; %#ok<AGROW>
            theta2DdotValues(kExp) = dx0(4);
            theta2ModelRangeValues(kExp) = theta2ModelRange;
            theta2MeasuredRangeValues(kExp) = theta2MeasuredRange;

            if makeZoomPlots && any(strcmp(stageName, stagesToPlot)) && any(strcmp(dataSetName, dataSetsToPlot))
                nexttile;
                plot(t(idxZoom), yMeasured(idxZoom, 2), 'LineWidth', 1.0, 'DisplayName', '\theta_2 measured');
                hold on;
                plot(t(idxZoom), yModel(idxZoom, 2), '--', 'LineWidth', 1.0, 'DisplayName', '\theta_2 model');
                grid on;
                title(sprintf('%s | ddot0 %.2f rad/s^2', expName, dx0(4)), 'Interpreter', 'none');
                xlabel('Time [s]');
                ylabel('\theta_2 [rad]');
                if kExp == 1
                    legend('Location', 'best');
                end
            end
        end

        fprintf('%s %s initial theta_2 acceleration summary:\n', stageName, dataSetName);
        fprintf('  mean(abs(theta2_ddot0))        = %.4g rad/s^2 = %.4g deg/s^2\n', ...
            mean(abs(theta2DdotValues), 'omitnan'), rad2deg(mean(abs(theta2DdotValues), 'omitnan')));
        fprintf('  max(abs(theta2_ddot0))         = %.4g rad/s^2 = %.4g deg/s^2\n', ...
            max(abs(theta2DdotValues)), rad2deg(max(abs(theta2DdotValues))));
        fprintf('  mean theta2 measured range %.1fs = %.4g rad = %.4g deg\n', ...
            restSecondsForMetrics, mean(theta2MeasuredRangeValues, 'omitnan'), rad2deg(mean(theta2MeasuredRangeValues, 'omitnan')));
        fprintf('  mean theta2 model range %.1fs    = %.4g rad = %.4g deg\n\n', ...
            restSecondsForMetrics, mean(theta2ModelRangeValues, 'omitnan'), rad2deg(mean(theta2ModelRangeValues, 'omitnan')));

        if makeZoomPlots && any(strcmp(stageName, stagesToPlot)) && any(strcmp(dataSetName, dataSetsToPlot))
            sgtitle(sprintf('%s %s: first %.1f s of \theta_2', stageName, dataSetName, zoomSeconds));
        end
    end
end

outputCsv = fullfile(scriptFolder, 'full_system_id_initial_theta2_oscillation_diagnostic_downline_60_40_theta2w15.csv');
writetable(allRows, outputCsv);

fprintf('Saved initial theta_2 diagnostic table:\n  %s\n\n', outputCsv);
fprintf('Interpretation guide:\n');
fprintf('  Large |Theta2Ddot0_rad_s2| means the model is not in static equilibrium at the measured initial state.\n');
fprintf('  If Theta2ModelRangeFirstRest is much larger than Theta2MeasuredRangeFirstRest,\n');
fprintf('  the early oscillation is caused by model acceleration during the initial rest segment.\n');
fprintf('  This can happen even when the initial velocity is exactly zero.\n');

%% ================================================================
% Local functions
%% ================================================================

function parameterValues = getParameterValues(model)
parameterValues = cell(length(model.Parameters), 1);
for k = 1:length(model.Parameters)
    parameterValues{k} = model.Parameters(k).Value;
end
end

function value = getParameterValueByName(model, name)
value = NaN;
for k = 1:length(model.Parameters)
    if strcmp(char(model.Parameters(k).Name), name)
        value = model.Parameters(k).Value;
        return;
    end
end
error('Could not find parameter %s.', name);
end

function zExp = getExperiment(zData, kExp)
try
    zExp = getexp(zData, kExp);
catch
    zExp = zData(:, :, kExp);
end
end

function modelOne = makeSingleExperimentModel(modelFile, order, parameterValues, x0, TsModel, modelTemplate)
initialStates = {x0(1); x0(2); x0(3); x0(4)};
modelOne = idnlgrey(modelFile, order, parameterValues, initialStates, TsModel);
for kPar = 1:length(modelTemplate.Parameters)
    modelOne.Parameters(kPar).Name = modelTemplate.Parameters(kPar).Name;
    modelOne.Parameters(kPar).Minimum = modelTemplate.Parameters(kPar).Minimum;
    modelOne.Parameters(kPar).Maximum = modelTemplate.Parameters(kPar).Maximum;
    modelOne.Parameters(kPar).Fixed = true;
end
modelOne.InputName = modelTemplate.InputName;
modelOne.InputUnit = modelTemplate.InputUnit;
modelOne.OutputName = modelTemplate.OutputName;
modelOne.OutputUnit = modelTemplate.OutputUnit;
modelOne.TimeUnit = modelTemplate.TimeUnit;
for kState = 1:length(modelTemplate.InitialStates)
    modelOne.InitialStates(kState).Name = modelTemplate.InitialStates(kState).Name;
    modelOne.InitialStates(kState).Unit = modelTemplate.InitialStates(kState).Unit;
end
modelOne = setinit(modelOne, 'Fixed', {true; true; true; true});
end

function y = getOutputArray(z)
try
    y = get(z, 'OutputData');
catch
    y = z.OutputData;
end
if iscell(y)
    y = y{1};
end
y = double(squeeze(y));
if size(y, 2) ~= 2 && size(y, 1) == 2
    y = y.';
end
if isvector(y) && size(y, 2) ~= 2
    y = y(:);
end
if size(y, 2) ~= 2
    error('Expected output data with 2 columns, but got size %s.', mat2str(size(y)));
end
end

function u = getInputArray(z)
try
    u = z.InputData;
catch
    try
        u = get(z, 'InputData');
    catch
        u = [];
    end
end
if iscell(u)
    u = u{1};
end
u = double(squeeze(u));
if isrow(u)
    u = u.';
end
end

function t = getTimeVector(z, N)
try
    Ts = z.Ts;
catch
    Ts = get(z, 'Ts');
end
if iscell(Ts)
    Ts = Ts{1};
end
try
    t = z.SamplingInstants;
catch
    try
        t = get(z, 'SamplingInstants');
    catch
        t = [];
    end
end
if iscell(t)
    t = t{1};
end
if isempty(t)
    t = (0:N-1).' * Ts;
else
    t = double(t(:));
    if length(t) ~= N
        t = (0:N-1).' * Ts;
    end
end
end
