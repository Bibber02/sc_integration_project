%% Kalman observer tuning
% User-facing script for tuning the observer process noise.

clc;

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = scriptFolder;
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(projectRoot);

startupFile = fullfile(projectRoot, 'startup.m');
if isfile(startupFile)
    run(startupFile);
else
    scip.setupPath;
end
projectPaths = scip.paths;

settings = kalman_default_settings(projectPaths);
if exist('kalmanTuningSettings', 'var') == 1 && isstruct(kalmanTuningSettings)
    settings = kalman_merge_settings(settings, kalmanTuningSettings);
end

if ~isfolder(settings.resultsFolder)
    mkdir(settings.resultsFolder);
end

fprintf('\n======================================================\n');
fprintf('Starting Kalman observer tuning\n');
fprintf('Reference runs: ');
fprintf('%d ', settings.idxReference);
fprintf('\nTuning runs:    ');
fprintf('%d ', settings.idxTune);
fprintf('\n======================================================\n\n');

sensorNoise = kalman_compute_sensor_covariance(settings.noiseDataFile);
R_kalman = sensorNoise.R;

fprintf('Measured sensor covariance R_kalman:\n');
disp(R_kalman);

observerModel = kalman_load_linear_observer_model(projectPaths, settings);
settings.runtimeSampleTime = observerModel.runtimeSampleTime;

allExperiments = kalman_load_full_system_experiments(projectPaths, settings);
referenceExperiments = kalman_select_experiments(allExperiments, settings.idxReference);
tuneExperiments = kalman_select_experiments(allExperiments, settings.idxTune);

if isempty(referenceExperiments)
    error('No reference experiments were loaded. Check idxReference and data settings.');
end
if isempty(tuneExperiments)
    error('No tuning experiments were loaded. Check idxTune and data settings.');
end

referenceID = kalman_identify_reference_model(projectPaths, settings, referenceExperiments);

referenceTune = kalman_simulate_reference_states(referenceID.model, tuneExperiments, observerModel, settings);
referenceHoldout = kalman_simulate_reference_states(referenceID.model, referenceExperiments, observerModel, settings);

P0_kalman = settings.P0_kalman;
x0_kalman = settings.x0_kalman;
defaultQRuntime = settings.defaultRuntimeQ;
defaultQOffline = defaultQRuntime;

normalization = kalman_build_normalization(referenceTune, R_kalman, settings);

fprintf('\nEvaluating current/default Kalman process noise...\n');
defaultTuneMetrics = kalman_evaluate_filter(defaultQOffline, R_kalman, P0_kalman, ...
    tuneExperiments, referenceTune, observerModel, settings, normalization);
defaultHoldoutMetrics = kalman_evaluate_filter(defaultQOffline, R_kalman, P0_kalman, ...
    referenceExperiments, referenceHoldout, observerModel, settings, normalization);

initialLogIntensities = settings.initialLogIntensities;
if isempty(initialLogIntensities)
    qAngle0 = max(defaultQRuntime(1, 1) / settings.runtimeSampleTime, eps);
    qRate0 = max(mean([defaultQRuntime(3, 3), defaultQRuntime(4, 4)]) / settings.runtimeSampleTime, eps);
    initialLogIntensities = log([qAngle0, qRate0]);
end

objective = @(logIntensities) kalman_tuning_objective(logIntensities, R_kalman, P0_kalman, ...
    tuneExperiments, referenceTune, observerModel, settings, normalization);

optimOptions = optimset( ...
    'Display', settings.optimizerDisplay, ...
    'MaxIter', settings.optimizerMaxIterations, ...
    'MaxFunEvals', settings.optimizerMaxFunctionEvaluations, ...
    'TolX', settings.optimizerTolerance, ...
    'TolFun', settings.optimizerTolerance);

fprintf('\nOptimizing process-noise intensities in log space...\n');
[logIntensitiesTuned, tunedScore] = fminsearch(objective, initialLogIntensities(:).', optimOptions);

processNoiseIntensities = exp(logIntensitiesTuned(:));
Q_tuned_offline = kalman_process_noise_matrix(processNoiseIntensities, settings.TsData);
Q_tuned_runtime = kalman_process_noise_matrix(processNoiseIntensities, settings.runtimeSampleTime);

tunedTuneMetrics = kalman_evaluate_filter(Q_tuned_offline, R_kalman, P0_kalman, ...
    tuneExperiments, referenceTune, observerModel, settings, normalization);
tunedHoldoutMetrics = kalman_evaluate_filter(Q_tuned_offline, R_kalman, P0_kalman, ...
    referenceExperiments, referenceHoldout, observerModel, settings, normalization);

if tunedTuneMetrics.meanScore < defaultTuneMetrics.meanScore
    Q_kalman = Q_tuned_runtime;
    selectedSource = 'tuned';
    selectedTuneMetrics = tunedTuneMetrics;
    selectedHoldoutMetrics = tunedHoldoutMetrics;
else
    Q_kalman = defaultQRuntime;
    selectedSource = 'default';
    selectedTuneMetrics = defaultTuneMetrics;
    selectedHoldoutMetrics = defaultHoldoutMetrics;
end

metrics = struct();
metrics.defaultTune = defaultTuneMetrics;
metrics.defaultHoldout = defaultHoldoutMetrics;
metrics.tunedTune = tunedTuneMetrics;
metrics.tunedHoldout = tunedHoldoutMetrics;
metrics.selectedTune = selectedTuneMetrics;
metrics.selectedHoldout = selectedHoldoutMetrics;
metrics.optimizedScore = tunedScore;
metrics.selectedSource = selectedSource;

kalmanSettings = settings;
referenceModelMetadata = referenceID.metadata;

resultFile = fullfile(settings.resultsFolder, 'kalman_tuning_result.mat');
save(resultFile, ...
    'Q_kalman', 'R_kalman', 'P0_kalman', 'x0_kalman', ...
    'processNoiseIntensities', 'Q_tuned_offline', 'Q_tuned_runtime', ...
    'defaultQOffline', 'defaultQRuntime', ...
    'metrics', 'sensorNoise', 'kalmanSettings', ...
    'referenceModelMetadata');

fprintf('\n======================================================\n');
fprintf('Kalman tuning complete\n');
fprintf('Selected process noise source: %s\n', selectedSource);
fprintf('Default tune score: %.6g\n', defaultTuneMetrics.meanScore);
fprintf('Tuned tune score:   %.6g\n', tunedTuneMetrics.meanScore);
fprintf('Saved result file:  %s\n', resultFile);
fprintf('======================================================\n\n');

disp('Selected runtime Q_kalman:');
disp(Q_kalman);
disp('Measured R_kalman:');
disp(R_kalman);

if settings.showPlots
    kalman_plot_tuning_metrics(metrics, settings);
end
