function Q = kalman_process_noise_matrix(processNoiseIntensities, sampleTime)
%KALMAN_PROCESS_NOISE_MATRIX Build diagonal discrete process covariance.

q = double(processNoiseIntensities(:));
if numel(q) ~= 2 || any(~isfinite(q)) || any(q <= 0)
    error('Process-noise intensities must be two positive finite values.');
end

Q = diag([q(1), q(1), q(2), q(2)]) * sampleTime;
end

