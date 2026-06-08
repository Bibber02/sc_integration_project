function experiments = kalman_load_full_system_experiments(projectPaths, settings)
%KALMAN_LOAD_FULL_SYSTEM_EXPERIMENTS Load PRBS/chirp full-system datasets.

experiments = struct( ...
    'name', {}, ...
    'signalType', {}, ...
    'runNumber', {}, ...
    'amplitude', {}, ...
    'file', {}, ...
    't', {}, ...
    'u', {}, ...
    'y', {}, ...
    'Ts', {});

dataFolder = projectPaths.fullSystemMeasurementData;
prbsFolder = fullfile(dataFolder, 'prbs');
chirpFolder = fullfile(dataFolder, 'chirp');

for kRun = 1:numel(settings.amplitudes)
    amplitude = settings.amplitudes(kRun);
    ampText = strrep(sprintf('%.2f', amplitude), '.', 'p');

    if settings.usePrbs
        filename = fullfile(prbsFolder, sprintf('fullsystem_prbs_A%s_run%02d.mat', ampText, kRun));
        experiments(end + 1) = loadOne(filename, 'prbs', kRun, amplitude, settings); %#ok<AGROW>
    end

    if settings.useChirp
        filename = fullfile(chirpFolder, sprintf('fullsystem_chirp_A%s_run%02d.mat', ampText, kRun));
        experiments(end + 1) = loadOne(filename, 'chirp', kRun, amplitude, settings); %#ok<AGROW>
    end
end
end

function experiment = loadOne(filename, signalType, runNumber, amplitude, settings)
if ~isfile(filename)
    error('Full-system data file not found: %s', filename);
end

S = load(filename, 'theta_1', 'theta_2', 'u_ts');
theta1 = kalman_timeseries_data(S.theta_1);
theta2 = kalman_timeseries_data(S.theta_2);
u = settings.inputSign * kalman_timeseries_data(S.u_ts);

if isa(S.theta_1, 'timeseries')
    t = double(squeeze(S.theta_1.Time(:)));
else
    t = (0:numel(theta1) - 1).' * settings.TsData;
end

n = min([numel(theta1), numel(theta2), numel(u), numel(t)]);
theta1 = theta1(1:n);
theta2 = theta2(1:n);
u = u(1:n);
t = t(1:n);

if n > 1
    Ts = median(diff(t));
else
    Ts = settings.TsData;
end

if abs(Ts - settings.TsData) > 1e-8
    warning('Experiment %s run %02d has Ts %.9g, expected %.9g.', ...
        signalType, runNumber, Ts, settings.TsData);
end

experiment = struct();
experiment.name = sprintf('%s_run_%02d', signalType, runNumber);
experiment.signalType = signalType;
experiment.runNumber = runNumber;
experiment.amplitude = amplitude;
experiment.file = filename;
experiment.t = t;
experiment.u = u;
experiment.y = [theta1, theta2];
experiment.Ts = Ts;
end

