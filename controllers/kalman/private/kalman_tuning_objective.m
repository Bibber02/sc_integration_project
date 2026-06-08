function score = kalman_tuning_objective(logIntensities, R, P0, experiments, references, observerModel, settings, normalization)
%KALMAN_TUNING_OBJECTIVE Objective minimized by fminsearch.

if any(~isfinite(logIntensities))
    score = Inf;
    return;
end

processNoiseIntensities = exp(logIntensities(:));
Q = kalman_process_noise_matrix(processNoiseIntensities, settings.TsData);
metrics = kalman_evaluate_filter(Q, R, P0, experiments, references, observerModel, settings, normalization);
score = metrics.meanScore;
end

