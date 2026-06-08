function sensorNoise = kalman_compute_sensor_covariance(noiseDataFile)
%KALMAN_COMPUTE_SENSOR_COVARIANCE Compute fixed R from static sensor data.

if ~isfile(noiseDataFile)
    error('Sensor noise data file not found: %s', noiseDataFile);
end

S = load(noiseDataFile, 'theta_1', 'theta_2');
if ~isfield(S, 'theta_1') || ~isfield(S, 'theta_2')
    error('Expected theta_1 and theta_2 in %s.', noiseDataFile);
end

theta1 = kalman_timeseries_data(S.theta_1);
theta2 = kalman_timeseries_data(S.theta_2);
n = min(numel(theta1), numel(theta2));
theta1 = theta1(1:n);
theta2 = theta2(1:n);

noise = [theta1 - mean(theta1), theta2 - mean(theta2)];
R = cov(noise, 1);
R = (R + R.') / 2;

if ~all(isfinite(R), 'all') || size(R, 1) ~= 2 || size(R, 2) ~= 2
    error('Computed measurement covariance R is invalid.');
end

[~, p] = chol(R);
if p ~= 0
    R = R + 1e-12 * eye(2);
    [~, p] = chol(R);
end
if p ~= 0
    error('Computed measurement covariance R is not positive definite.');
end

sensorNoise = struct();
sensorNoise.file = noiseDataFile;
sensorNoise.R = R;
sensorNoise.mean = [mean(theta1); mean(theta2)];
sensorNoise.std = sqrt(diag(R));
sensorNoise.numSamples = n;
sensorNoise.correlation = R(1, 2) / sqrt(R(1, 1) * R(2, 2));
end

