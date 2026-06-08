clear;
clc;
close all;

%% ================================================================
% Full-system black-box identification using PRBS and chirp data
%
% This script is the black-box counterpart to the grey-box full-system
% identification script. It does not use the nonlinear equations of motion
% or the passive-link parameter file. Instead, it estimates a discrete-time
% linear MIMO state-space model directly from measured input/output data.
%
% Data expectation per .mat file:
%   theta_1   timeseries with measured rod-1 angle [rad]
%   theta_2   timeseries with measured rod-2 angle [rad]
%   u_ts      timeseries with input voltage/control signal [V]
%
% The script:
%   1) loads PRBS and chirp experiments,
%   2) builds multi-experiment iddata objects,
%   3) optionally removes constant offsets for deviation-variable ID,
%   4) estimates state-space models for several candidate orders,
%   5) selects the order with the best validation fit,
%   6) saves the best black-box model and comparison results.
% ================================================================

%% ================================================================
% Configuration
% ================================================================

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = scriptFolder;
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(projectRoot);
scip.setupPath;
projectPaths = scip.paths;

% Folder layout. Change these paths if your project structure is different.
dataFolder   = projectPaths.fullSystemMeasurementData;
prbsFolder   = fullfile(dataFolder, 'prbs');
chirpFolder  = fullfile(dataFolder, 'chirp');
outputFolder = scriptFolder;

% Data settings.
Ts = 0.01;
inputSign = -1;             % Match the sign convention used in the grey-box script.
amplitudes = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];
idxEst = [1 3 5 7 9];       % Odd runs are used for identification.
idxVal = [2 4 6 8 10];      % Even runs are used for validation.
usePRBS = true;
useChirp = true;

% Preprocessing.
% For a black-box linear model it is usually better to identify the dynamics
% in deviation variables around the operating point. Therefore, removing the
% constant mean is enabled by default. Turn this off if you explicitly want
% to fit the absolute measured angles instead.
removeConstantOffsets = true;

% Black-box model settings.
modelOrders = 2:8;                  % Try 4 first mentally, but validate several orders.
estimationFocus = 'simulation';     % 'simulation' is usually better for open-loop simulation fit.
enforceStability = true;
forceNoDirectFeedthrough = true;    % Physical plant has no instantaneous voltage-to-angle feedthrough.
maxIterations = 150;

% Output and plot settings.
showRawDataPlots = true;
showComparePlots = true;
showBodePlot = true;
showPoleZeroPlot = true;
figurePosition = [80 80 1050 650];

consoleLogFile = fullfile(outputFolder, 'full_system_blackbox_console_output.txt');
resultMatFile = fullfile(outputFolder, 'full_system_blackbox_state_space_result.mat');
resultCsvFile = fullfile(outputFolder, 'full_system_blackbox_state_space_order_results.csv');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

diary(consoleLogFile);

fprintf('\n======================================================\n');
fprintf('Starting full-system BLACK-BOX state-space identification\n');
fprintf('No grey-box equations and no passive-link parameters are used\n');
fprintf('======================================================\n\n');

fprintf('Identification runs: 1, 3, 5, 7, 9 for selected signal types\n');
fprintf('Validation runs:     2, 4, 6, 8, 10 for selected signal types\n');
fprintf('removeConstantOffsets = %d\n\n', removeConstantOffsets);

%% ================================================================
% Load data and build multi-experiment iddata objects
% ================================================================

zEstRaw = [];
zValRaw = [];

nEst = 0;
nVal = 0;

experimentNamesEst = {};
experimentNamesVal = {};

rawPlotData = struct('name', {}, 't', {}, 'theta1', {}, 'theta2', {}, 'u', {}, 'run', {}, 'type', {});
rawPlotCounter = 0;

for k = 1:length(amplitudes)

    ampText = strrep(sprintf('%.2f', amplitudes(k)), '.', 'p');

    if usePRBS
        prbsFile = fullfile(prbsFolder, sprintf('fullsystem_prbs_A%s_run%02d.mat', ampText, k));

        if ~isfile(prbsFile)
            error('Could not find PRBS file: %s', prbsFile);
        end

        load(prbsFile, 'theta_1', 'theta_2', 'u_ts');

        theta1 = double(squeeze(theta_1.Data(:)));
        theta2 = double(squeeze(theta_2.Data(:)));
        u = inputSign * double(squeeze(u_ts.Data(:)));
        t = double(squeeze(theta_1.Time(:)));

        nMin = min([length(theta1), length(theta2), length(u)]);
        theta1 = theta1(1:nMin);
        theta2 = theta2(1:nMin);
        u = u(1:nMin);
        if length(t) >= nMin
            t = t(1:nMin);
        else
            t = (0:nMin-1).' * Ts;
        end

        y = [theta1, theta2];
        z = iddata(y, u, Ts);
        z.Name = sprintf('prbs_run_%02d', k);
        z.InputName = {'u'};
        z.InputUnit = {'V'};
        z.OutputName = {'theta_1', 'theta_2'};
        z.OutputUnit = {'rad', 'rad'};
        z.TimeUnit = 's';

        if ismember(k, idxEst)
            nEst = nEst + 1;
            if isempty(zEstRaw)
                zEstRaw = z;
            else
                zEstRaw = merge(zEstRaw, z);
            end
            experimentNamesEst{nEst, 1} = z.Name; %#ok<SAGROW>
        elseif ismember(k, idxVal)
            nVal = nVal + 1;
            if isempty(zValRaw)
                zValRaw = z;
            else
                zValRaw = merge(zValRaw, z);
            end
            experimentNamesVal{nVal, 1} = z.Name; %#ok<SAGROW>
        end

        rawPlotCounter = rawPlotCounter + 1;
        rawPlotData(rawPlotCounter).name = z.Name;
        rawPlotData(rawPlotCounter).t = t;
        rawPlotData(rawPlotCounter).theta1 = theta1;
        rawPlotData(rawPlotCounter).theta2 = theta2;
        rawPlotData(rawPlotCounter).u = u;
        rawPlotData(rawPlotCounter).run = k;
        rawPlotData(rawPlotCounter).type = 'PRBS';
    end

    if useChirp
        chirpFile = fullfile(chirpFolder, sprintf('fullsystem_chirp_A%s_run%02d.mat', ampText, k));

        if ~isfile(chirpFile)
            error('Could not find chirp file: %s', chirpFile);
        end

        load(chirpFile, 'theta_1', 'theta_2', 'u_ts');

        theta1 = double(squeeze(theta_1.Data(:)));
        theta2 = double(squeeze(theta_2.Data(:)));
        u = inputSign * double(squeeze(u_ts.Data(:)));
        t = double(squeeze(theta_1.Time(:)));

        nMin = min([length(theta1), length(theta2), length(u)]);
        theta1 = theta1(1:nMin);
        theta2 = theta2(1:nMin);
        u = u(1:nMin);
        if length(t) >= nMin
            t = t(1:nMin);
        else
            t = (0:nMin-1).' * Ts;
        end

        y = [theta1, theta2];
        z = iddata(y, u, Ts);
        z.Name = sprintf('chirp_run_%02d', k);
        z.InputName = {'u'};
        z.InputUnit = {'V'};
        z.OutputName = {'theta_1', 'theta_2'};
        z.OutputUnit = {'rad', 'rad'};
        z.TimeUnit = 's';

        if ismember(k, idxEst)
            nEst = nEst + 1;
            if isempty(zEstRaw)
                zEstRaw = z;
            else
                zEstRaw = merge(zEstRaw, z);
            end
            experimentNamesEst{nEst, 1} = z.Name; %#ok<SAGROW>
        elseif ismember(k, idxVal)
            nVal = nVal + 1;
            if isempty(zValRaw)
                zValRaw = z;
            else
                zValRaw = merge(zValRaw, z);
            end
            experimentNamesVal{nVal, 1} = z.Name; %#ok<SAGROW>
        end

        rawPlotCounter = rawPlotCounter + 1;
        rawPlotData(rawPlotCounter).name = z.Name;
        rawPlotData(rawPlotCounter).t = t;
        rawPlotData(rawPlotCounter).theta1 = theta1;
        rawPlotData(rawPlotCounter).theta2 = theta2;
        rawPlotData(rawPlotCounter).u = u;
        rawPlotData(rawPlotCounter).run = k;
        rawPlotData(rawPlotCounter).type = 'Chirp';
    end
end

fprintf('Number of identification experiments: %d\n', nEst);
fprintf('Number of validation experiments:     %d\n\n', nVal);

if isempty(zEstRaw) || isempty(zValRaw)
    error('No estimation or validation data was loaded. Check usePRBS/useChirp and folder paths.');
end

%% ================================================================
% Optional raw-data plots
% ================================================================

if showRawDataPlots
    figure('Name', 'Black-box ID: raw full-system measurements', 'Position', figurePosition);
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    hold on;
    for i = 1:numel(rawPlotData)
        plot(rawPlotData(i).t, rawPlotData(i).theta1, 'DisplayName', rawPlotData(i).name);
    end
    title('\theta_1 measurements');
    xlabel('Time [s]');
    ylabel('\theta_1 [rad]');
    grid on;

    nexttile;
    hold on;
    for i = 1:numel(rawPlotData)
        plot(rawPlotData(i).t, rawPlotData(i).theta2, 'DisplayName', rawPlotData(i).name);
    end
    title('\theta_2 measurements');
    xlabel('Time [s]');
    ylabel('\theta_2 [rad]');
    grid on;

    nexttile;
    hold on;
    for i = 1:numel(rawPlotData)
        plot(rawPlotData(i).t, rawPlotData(i).u, 'DisplayName', rawPlotData(i).name);
    end
    title('Input signal');
    xlabel('Time [s]');
    ylabel('u [V]');
    grid on;
end

%% ================================================================
% Preprocess data for black-box identification
% ================================================================

zEst = zEstRaw;
zVal = zValRaw;

if removeConstantOffsets
    zEst = detrend(zEstRaw, 0);
    zVal = detrend(zValRaw, 0);
    fprintf('Constant offsets removed from each experiment using detrend(data,0).\n');
else
    fprintf('Constant offsets were NOT removed. The model is fitted to absolute signal levels.\n');
end

%% ================================================================
% Black-box state-space model estimation
% ================================================================

compareOpt = compareOptions;
compareOpt.InitialCondition = 'estimate';

optN4 = n4sidOptions;
optN4.Display = 'off';
optN4.Focus = estimationFocus;
optN4.InitialState = 'estimate';
optN4.EnforceStability = enforceStability;

optSS = ssestOptions;
optSS.Display = 'off';
optSS.Focus = estimationFocus;
optSS.InitialState = 'estimate';
optSS.EnforceStability = enforceStability;
optSS.SearchOptions.MaxIterations = maxIterations;

nOrders = numel(modelOrders);
allModels = cell(nOrders, 1);
allN4Models = cell(nOrders, 1);

orderResult = zeros(nOrders, 5);   % [order, meanN4Est, meanN4Val, meanSSEst, meanSSVal]

bestValidationFit = -Inf;
bestOrder = NaN;
modelBest = [];
modelBestN4SID = [];

for io = 1:nOrders

    nx = modelOrders(io);

    fprintf('\n======================================================\n');
    fprintf('Estimating black-box state-space model of order %d\n', nx);
    fprintf('======================================================\n\n');

    orderResult(io, :) = [nx, NaN, NaN, NaN, NaN];

    try
        modelN4 = n4sid(zEst, nx, optN4);
        modelN4.Name = sprintf('N4SID order %d', nx);
        modelN4.InputName = {'u'};
        modelN4.InputUnit = {'V'};
        modelN4.OutputName = {'theta_1', 'theta_2'};
        modelN4.OutputUnit = {'rad', 'rad'};
        modelN4.TimeUnit = 's';

        if forceNoDirectFeedthrough
            try
                modelN4.D(:) = 0;
                modelN4.Structure.D.Free(:) = false;
            catch ME_D
                warning('Could not fix D matrix for N4SID order %d: %s', nx, ME_D.message);
            end
        end

        allN4Models{io} = modelN4;

        [~, fitN4Est] = compare(zEst, modelN4, compareOpt);
        [~, fitN4Val] = compare(zVal, modelN4, compareOpt);

        fitN4EstValues = [];
        if iscell(fitN4Est)
            for kk = 1:numel(fitN4Est)
                fitN4EstValues = [fitN4EstValues; fitN4Est{kk}(:)]; %#ok<AGROW>
            end
        else
            fitN4EstValues = fitN4Est(:);
        end
        fitN4EstValues = fitN4EstValues(isfinite(fitN4EstValues));

        fitN4ValValues = [];
        if iscell(fitN4Val)
            for kk = 1:numel(fitN4Val)
                fitN4ValValues = [fitN4ValValues; fitN4Val{kk}(:)]; %#ok<AGROW>
            end
        else
            fitN4ValValues = fitN4Val(:);
        end
        fitN4ValValues = fitN4ValValues(isfinite(fitN4ValValues));

        meanN4Est = mean(fitN4EstValues);
        meanN4Val = mean(fitN4ValValues);

        fprintf('N4SID order %d mean estimation fit: %.2f %%\n', nx, meanN4Est);
        fprintf('N4SID order %d mean validation fit: %.2f %%\n', nx, meanN4Val);

        try
            modelSS = ssest(zEst, modelN4, optSS);
        catch ME_SS_INIT
            warning('ssest refinement from N4SID failed for order %d: %s', nx, ME_SS_INIT.message);
            warning('Trying direct ssest(data,order,opt) for order %d.', nx);
            modelSS = ssest(zEst, nx, optSS);
        end

        modelSS.Name = sprintf('SSEST refined order %d', nx);
        modelSS.InputName = {'u'};
        modelSS.InputUnit = {'V'};
        modelSS.OutputName = {'theta_1', 'theta_2'};
        modelSS.OutputUnit = {'rad', 'rad'};
        modelSS.TimeUnit = 's';

        if forceNoDirectFeedthrough
            try
                modelSS.D(:) = 0;
                modelSS.Structure.D.Free(:) = false;
            catch ME_D2
                warning('Could not keep D matrix fixed for SSEST order %d: %s', nx, ME_D2.message);
            end
        end

        allModels{io} = modelSS;

        [~, fitSSEst] = compare(zEst, modelSS, compareOpt);
        [~, fitSSVal] = compare(zVal, modelSS, compareOpt);

        fitSSEstValues = [];
        if iscell(fitSSEst)
            for kk = 1:numel(fitSSEst)
                fitSSEstValues = [fitSSEstValues; fitSSEst{kk}(:)]; %#ok<AGROW>
            end
        else
            fitSSEstValues = fitSSEst(:);
        end
        fitSSEstValues = fitSSEstValues(isfinite(fitSSEstValues));

        fitSSValValues = [];
        if iscell(fitSSVal)
            for kk = 1:numel(fitSSVal)
                fitSSValValues = [fitSSValValues; fitSSVal{kk}(:)]; %#ok<AGROW>
            end
        else
            fitSSValValues = fitSSVal(:);
        end
        fitSSValValues = fitSSValValues(isfinite(fitSSValValues));

        meanSSEst = mean(fitSSEstValues);
        meanSSVal = mean(fitSSValValues);

        orderResult(io, :) = [nx, meanN4Est, meanN4Val, meanSSEst, meanSSVal];

        fprintf('SSEST order %d mean estimation fit: %.2f %%\n', nx, meanSSEst);
        fprintf('SSEST order %d mean validation fit: %.2f %%\n', nx, meanSSVal);

        if meanSSVal > bestValidationFit
            bestValidationFit = meanSSVal;
            bestOrder = nx;
            modelBest = modelSS;
            modelBestN4SID = modelN4;
        end

    catch ME
        warning('Order %d failed: %s', nx, ME.message);
    end
end

orderResultsTable = array2table(orderResult, ...
    'VariableNames', {'Order', 'MeanN4SIDEstFit', 'MeanN4SIDValFit', 'MeanSSESTEstFit', 'MeanSSESTValFit'});

disp(orderResultsTable);
writetable(orderResultsTable, resultCsvFile);

if isempty(modelBest)
    error('All black-box model orders failed. Check data, options, and System Identification Toolbox availability.');
end

fprintf('\n======================================================\n');
fprintf('Selected black-box model order: %d\n', bestOrder);
fprintf('Selected validation fit: %.2f %%\n', bestValidationFit);
fprintf('======================================================\n\n');

%% ================================================================
% Final comparisons and plots
% ================================================================

[yEstBest, fitEstBest] = compare(zEst, modelBest, compareOpt);
[yValBest, fitValBest] = compare(zVal, modelBest, compareOpt);

fitEstValues = [];
if iscell(fitEstBest)
    for kk = 1:numel(fitEstBest)
        fitEstValues = [fitEstValues; fitEstBest{kk}(:)]; %#ok<AGROW>
    end
else
    fitEstValues = fitEstBest(:);
end
fitEstValues = fitEstValues(isfinite(fitEstValues));

fitValValues = [];
if iscell(fitValBest)
    for kk = 1:numel(fitValBest)
        fitValValues = [fitValValues; fitValBest{kk}(:)]; %#ok<AGROW>
    end
else
    fitValValues = fitValBest(:);
end
fitValValues = fitValValues(isfinite(fitValValues));

fprintf('Best model final mean estimation fit: %.2f %%\n', mean(fitEstValues));
fprintf('Best model final mean validation fit: %.2f %%\n\n', mean(fitValValues));

fprintf('Final identification fit, per experiment and output:\n');
disp(fitEstBest);

fprintf('Final validation fit, per experiment and output:\n');
disp(fitValBest);

if showComparePlots
    figure('Name', 'Black-box final identification fit', 'Position', figurePosition);
    compare(zEst, modelBest, compareOpt);
    title(sprintf('Black-box state-space identification fit, order %d', bestOrder));

    figure('Name', 'Black-box final validation fit', 'Position', figurePosition + [30 30 0 0]);
    compare(zVal, modelBest, compareOpt);
    title(sprintf('Black-box state-space validation fit, order %d', bestOrder));
end

if showBodePlot
    figure('Name', 'Black-box model Bode plot', 'Position', figurePosition + [60 60 0 0]);
    bode(modelBest);
    grid on;
    title(sprintf('Black-box state-space Bode plot, order %d', bestOrder));
end

if showPoleZeroPlot
    figure('Name', 'Black-box model pole-zero map', 'Position', figurePosition + [90 90 0 0]);
    pzmap(modelBest);
    grid on;
    title(sprintf('Black-box state-space pole-zero map, order %d', bestOrder));
end

%% ================================================================
% Save results for Simulink/control use
% ================================================================

% modelBest is an identified state-space model. sysdBest is the plain LTI
% state-space model that can be used in standard Control System Toolbox blocks.
sysdBest = ss(modelBest);

save(resultMatFile, ...
    'modelBest', 'modelBestN4SID', 'sysdBest', 'allModels', 'allN4Models', ...
    'zEst', 'zVal', 'zEstRaw', 'zValRaw', 'yEstBest', 'yValBest', ...
    'fitEstBest', 'fitValBest', 'orderResultsTable', 'bestOrder', ...
    'bestValidationFit', 'modelOrders', 'Ts', 'inputSign', 'amplitudes', ...
    'idxEst', 'idxVal', 'removeConstantOffsets', 'experimentNamesEst', ...
    'experimentNamesVal');

fprintf('\nSaved result file: %s\n', resultMatFile);
fprintf('Saved order-result table: %s\n', resultCsvFile);
fprintf('\nUse sysdBest as the discrete-time LTI model in Simulink/control design.\n');

if removeConstantOffsets
    fprintf('Note: sysdBest is a deviation-variable model because constant offsets were removed.\n');
end

diary off;
