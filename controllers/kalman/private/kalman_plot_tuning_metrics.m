function kalman_plot_tuning_metrics(metrics, settings)
%KALMAN_PLOT_TUNING_METRICS Plot default-vs-tuned scores.

figure('Name', 'Kalman tuning metrics', 'Position', settings.figurePosition);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(categorical({'Default tune', 'Tuned tune', 'Default holdout', 'Tuned holdout'}), ...
    [metrics.defaultTune.meanScore, metrics.tunedTune.meanScore, ...
    metrics.defaultHoldout.meanScore, metrics.tunedHoldout.meanScore]);
ylabel('Normalized score');
grid on;

nexttile;
bar(categorical({'Default tune', 'Tuned tune', 'Default holdout', 'Tuned holdout'}), ...
    [metrics.defaultTune.meanAngleRmse, metrics.tunedTune.meanAngleRmse, ...
    metrics.defaultHoldout.meanAngleRmse, metrics.tunedHoldout.meanAngleRmse]);
ylabel('Mean angle RMSE [rad]');
grid on;
end

