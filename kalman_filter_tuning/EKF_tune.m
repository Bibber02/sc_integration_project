clear;
clc;
close all;

%% Configuration
scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptFolder);
addpath(fullfile(projectRoot, 'model'));

Ts = 0.01;
validationRuns = [2 4 6 8 10];
amplitudes = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];
usePrbs = true;
useChirp = true;
inputSign = -1;

P0 = 10 * eye(4);
initialQdiag = [0.1 0.1 1 0.1];
logQLower = [-10 -10 -8 -8];
logQUpper = [1 1 3 3];
nRandomStarts = 12;
randomSeed = 23;
outOfRangePenaltyWeight = 1e6;
plotRuns = [2 10];

%% Load model parameters and validation data
p = load_parameters();
experiments = loadValidationExperiments(projectRoot, Ts, amplitudes, ...
    validationRuns, usePrbs, useChirp, inputSign);

%% Fixed measurement-noise covariance
noiseFile = fullfile(projectRoot, ...
    'system_identification', 'sensor_noise_measurement', 'noise_data.mat');
noiseData = load(noiseFile, 'theta_1', 'theta_2');
noiseTheta1 = columnData(noiseData.theta_1);
noiseTheta2 = columnData(noiseData.theta_2);
nNoise = min(numel(noiseTheta1), numel(noiseTheta2));
noiseTheta1 = noiseTheta1(1:nNoise);
noiseTheta2 = noiseTheta2(1:nNoise);

noise = [noiseTheta1 - mean(noiseTheta1), noiseTheta2 - mean(noiseTheta2)];
R_ekf = cov(noise, 1);
R_ekf = (R_ekf + R_ekf.') / 2;

%% Q parameterization and EKF functions
stateTransitionFcn = @(x, u) rk4Transition(x, u, p, Ts);
measurementFcn = @measureAngles;
costFcn = @(logQdiag) ekfValidationCost(logQdiag, experiments, ...
    stateTransitionFcn, measurementFcn, R_ekf, P0, ...
    logQLower, logQUpper, outOfRangePenaltyWeight);

%% Optimize Q in log space with multiple starts
initialLogQdiag = log10(initialQdiag);
options = optimset( ...
    'Display', 'iter', ...
    'MaxIter', 80, ...
    'MaxFunEvals', 200, ...
    'TolX', 1e-3, ...
    'TolFun', 1e-3);

rng(randomSeed);
randomStartLogQdiag = logQLower + rand(nRandomStarts, 4) .* ...
    (logQUpper - logQLower);
startLogQdiag = [initialLogQdiag; randomStartLogQdiag];
nStarts = size(startLogQdiag, 1);

startResults = table((1:nStarts).', NaN(nStarts, 1), NaN(nStarts, 1), ...
    NaN(nStarts, 1), NaN(nStarts, 1), NaN(nStarts, 1), ...
    'VariableNames', {'Start', 'FinalCost', 'log_q_theta1', ...
    'log_q_theta2', 'log_q_omega1', 'log_q_omega2'});

bestLogQdiag = [];
finalValidationCost = inf;

for kStart = 1:nStarts
    fprintf('\nMulti-start EKF Q optimization %d of %d\n', kStart, nStarts);
    candidateLogQdiag = fminsearch(costFcn, ...
        startLogQdiag(kStart, :), options);

    candidateLogQdiag = clamp(candidateLogQdiag(:).', logQLower, logQUpper);
    candidateCost = costFcn(candidateLogQdiag);

    startResults.FinalCost(kStart) = candidateCost;
    startResults{kStart, 3:6} = candidateLogQdiag;

    if candidateCost < finalValidationCost
        finalValidationCost = candidateCost;
        bestLogQdiag = candidateLogQdiag;
    end
end

bestQdiag = 10 .^ bestLogQdiag(:).';

best_q_theta1 = bestQdiag(1);
best_q_theta2 = bestQdiag(2);
best_q_omega1 = bestQdiag(3);
best_q_omega2 = bestQdiag(4);
Q_ekf = diag(bestQdiag);

fprintf('\nBest EKF process-noise values:\n');
fprintf('  q_theta1 = %.12g\n', best_q_theta1);
fprintf('  q_theta2 = %.12g\n', best_q_theta2);
fprintf('  q_omega1 = %.12g\n', best_q_omega1);
fprintf('  q_omega2 = %.12g\n', best_q_omega2);
fprintf('\nQ_ekf:\n');
disp(Q_ekf);
fprintf('R_ekf:\n');
disp(R_ekf);
fprintf('Final average validation cost: %.12g\n', finalValidationCost);
fprintf('\nMulti-start summary:\n');
disp(sortrows(startResults, 'FinalCost'));

%% Final EKF run and plots
finalResults = cell(size(experiments));
for kExp = 1:numel(experiments)
    finalResults{kExp} = runEkfExperiment(experiments(kExp), Q_ekf, R_ekf, ...
        P0, stateTransitionFcn, measurementFcn);
end

for kExp = 1:numel(experiments)
    expData = experiments(kExp);
    if ismember(expData.runNumber, plotRuns)
        plotEkfExperiment(expData, finalResults{kExp});
    end
end

%% Local functions
function experiments = loadValidationExperiments(projectRoot, Ts, amplitudes, ...
    validationRuns, usePrbs, useChirp, inputSign)

dataFolder = fullfile(projectRoot, 'system_identification', ...
    'full_system', 'measurement_data');
prbsFolder = fullfile(dataFolder, 'prbs');
chirpFolder = fullfile(dataFolder, 'chirp');

experiments = struct( ...
    'name', {}, ...
    'signalType', {}, ...
    'runNumber', {}, ...
    'amplitude', {}, ...
    't', {}, ...
    'u', {}, ...
    'y', {});

for kRun = validationRuns
    amplitude = amplitudes(kRun);
    ampText = strrep(sprintf('%.2f', amplitude), '.', 'p');

    if usePrbs
        filename = fullfile(prbsFolder, ...
            sprintf('fullsystem_prbs_A%s_run%02d.mat', ampText, kRun));
        experiments(end + 1) = loadOneExperiment(filename, 'prbs', ...
            kRun, amplitude, Ts, inputSign); %#ok<AGROW>
    end

    if useChirp
        filename = fullfile(chirpFolder, ...
            sprintf('fullsystem_chirp_A%s_run%02d.mat', ampText, kRun));
        experiments(end + 1) = loadOneExperiment(filename, 'chirp', ...
            kRun, amplitude, Ts, inputSign); %#ok<AGROW>
    end
end
end

function experiment = loadOneExperiment(filename, signalType, runNumber, ...
    amplitude, Ts, inputSign)

S = load(filename, 'theta_1', 'theta_2', 'u_ts');
theta1 = columnData(S.theta_1);
theta2 = columnData(S.theta_2);
u = inputSign * columnData(S.u_ts);

if isa(S.theta_1, 'timeseries')
    t = double(squeeze(S.theta_1.Time(:)));
else
    t = (0:numel(theta1) - 1).' * Ts;
end

n = min([numel(theta1), numel(theta2), numel(u), numel(t)]);

experiment = struct();
experiment.name = sprintf('%s_run_%02d', signalType, runNumber);
experiment.signalType = signalType;
experiment.runNumber = runNumber;
experiment.amplitude = amplitude;
experiment.t = t(1:n);
experiment.u = u(1:n);
experiment.y = [theta1(1:n), theta2(1:n)];
end

function data = columnData(value)
if isa(value, 'timeseries')
    data = value.Data;
else
    data = value;
end

data = double(squeeze(data));
data = data(:);
end

function xNext = rk4Transition(x, u, p, Ts)
k1 = nonlinearPlant(x, u, p);
k2 = nonlinearPlant(x + 0.5 * Ts * k1, u, p);
k3 = nonlinearPlant(x + 0.5 * Ts * k2, u, p);
k4 = nonlinearPlant(x + Ts * k3, u, p);
xNext = x + (Ts / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
end

function y = measureAngles(x)
y = x(1:2);
end

function cost = ekfValidationCost(logQdiag, experiments, stateTransitionFcn, ...
    measurementFcn, R, P0, logQLower, logQUpper, penaltyWeight)

logQdiag = logQdiag(:).';
if any(~isfinite(logQdiag))
    cost = realmax;
    return;
end

lowerViolation = max(logQLower - logQdiag, 0);
upperViolation = max(logQdiag - logQUpper, 0);
penalty = penaltyWeight * sum(lowerViolation.^2 + upperViolation.^2);

logQdiag = clamp(logQdiag, logQLower, logQUpper);
qDiag = 10 .^ logQdiag;
Q = diag(qDiag);
variance = diag(R).';
costSum = 0;
nSamples = 0;

for kExp = 1:numel(experiments)
    result = runEkfExperiment(experiments(kExp), Q, R, P0, ...
        stateTransitionFcn, measurementFcn);
    innovation = result.innovation(2:end, :);
    normalizedInnovation = innovation(:, 1).^2 / variance(1) + ...
        innovation(:, 2).^2 / variance(2);
    costSum = costSum + sum(normalizedInnovation);
    nSamples = nSamples + numel(normalizedInnovation);
end

cost = costSum / nSamples + penalty;
end

function value = clamp(value, lower, upper)
value = min(max(value, lower), upper);
end

function result = runEkfExperiment(expData, Q, R, P0, stateTransitionFcn, ...
    measurementFcn)

y = expData.y;
u = expData.u;
n = size(y, 1);

x0 = [y(1, 1); y(1, 2); 0; 0];
ekf = extendedKalmanFilter(stateTransitionFcn, measurementFcn, x0);
ekf.ProcessNoise = Q;
ekf.MeasurementNoise = R;
ekf.StateCovariance = P0;

xPred = zeros(n, 4);
xEst = zeros(n, 4);
yPred = zeros(n, 2);
yEst = zeros(n, 2);
innovation = NaN(n, 2);

xPred(1, :) = x0.';
xEst(1, :) = x0.';
yPred(1, :) = measurementFcn(x0).';
yEst(1, :) = yPred(1, :);

for k = 2:n
    xPredK = predict(ekf, u(k - 1));
    yPredK = measurementFcn(xPredK);
    yK = y(k, :).';
    innovationK = yK - yPredK;

    correct(ekf, yK);
    xEstK = ekf.State;
    yEstK = measurementFcn(xEstK);

    xPred(k, :) = xPredK(:).';
    xEst(k, :) = xEstK(:).';
    yPred(k, :) = yPredK(:).';
    yEst(k, :) = yEstK(:).';
    innovation(k, :) = innovationK(:).';
end

result = struct();
result.xPred = xPred;
result.xEst = xEst;
result.yPred = yPred;
result.yEst = yEst;
result.innovation = innovation;
end

function plotEkfExperiment(expData, result)
t = expData.t;

figure('Name', sprintf('EKF output comparison: %s', expData.name));
tiledlayout(2, 1);

nexttile;
plot(t, expData.y(:, 1), 'k', ...
    t, result.yPred(:, 1), 'r--', ...
    t, result.yEst(:, 1), 'b');
grid on;
ylabel('\theta_1 [rad]');
title(sprintf('%s, A = %.2f', expData.name, expData.amplitude), ...
    'Interpreter', 'none');
legend('measured', 'predicted before correction', 'estimated after correction');

nexttile;
plot(t, expData.y(:, 2), 'k', ...
    t, result.yPred(:, 2), 'r--', ...
    t, result.yEst(:, 2), 'b');
grid on;
xlabel('Time [s]');
ylabel('\theta_2 [rad]');

figure('Name', sprintf('EKF estimated rates: %s', expData.name));
plot(t, result.xEst(:, 3), 'b', ...
    t, result.xEst(:, 4), 'r');
grid on;
xlabel('Time [s]');
ylabel('Estimated angular rate [rad/s]');
title(sprintf('Estimated rates: %s', expData.name), 'Interpreter', 'none');
legend('\omega_1 estimate', '\omega_2 estimate');
end
