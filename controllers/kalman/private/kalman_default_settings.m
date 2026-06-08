function settings = kalman_default_settings(projectPaths)
%KALMAN_DEFAULT_SETTINGS Defaults for Kalman observer tuning.

settings = struct();

settings.noiseDataFile = fullfile(projectPaths.sensorNoiseMeasurement, 'noise_data.mat');
settings.resultsFolder = fullfile(projectPaths.kalman, 'results');
settings.referenceResultFile = fullfile(settings.resultsFolder, 'kalman_reference_model.mat');

settings.TsData = 0.01;
settings.runtimeSampleTime = 0.001;
settings.inputSign = -1;

settings.amplitudes = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];
settings.idxReference = [2 4 6 8 10];
settings.idxTune = [1 3 5 7 9];
settings.usePrbs = true;
settings.useChirp = true;

settings.forceReferenceIdentification = false;
settings.referenceStage1MaxIterations = 15;
settings.referenceStage2MaxIterations = 25;
settings.referenceDisplay = 'Full';
settings.estimateReferenceCovariance = false;
settings.referenceOutputWeight = diag([1, 8]);

settings.optimizerMaxIterations = 40;
settings.optimizerMaxFunctionEvaluations = 80;
settings.optimizerTolerance = 1e-3;
settings.optimizerDisplay = 'iter';
settings.initialLogIntensities = [];

settings.defaultRuntimeQ = diag([0.1, 0.1, 1, 0.1]);
settings.P0_kalman = 10 * eye(4);
settings.x0_kalman = zeros(4, 1);

settings.initialEstimateMode = 'measured_angles';
settings.warmupSeconds = 2.0;
settings.minRateScale = 1e-3;

settings.showPlots = usejava('desktop');
settings.figurePosition = [80 80 1100 650];
end
