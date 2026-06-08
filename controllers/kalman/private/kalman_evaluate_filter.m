function metrics = kalman_evaluate_filter(Q, R, P0, experiments, references, observerModel, settings, normalization)
%KALMAN_EVALUATE_FILTER Score a Kalman filter against reference states.

nExp = numel(experiments);
names = strings(nExp, 1);
score = zeros(nExp, 1);
theta1Rmse = zeros(nExp, 1);
theta2Rmse = zeros(nExp, 1);
rate1Rmse = zeros(nExp, 1);
rate2Rmse = zeros(nExp, 1);
innovationRmse1 = zeros(nExp, 1);
innovationRmse2 = zeros(nExp, 1);

for i = 1:nExp
    expData = experiments(i);
    ref = references(i);
    names(i) = string(expData.name);

    yDev = expData.y - observerModel.y0(:).';
    uDev = expData.u - observerModel.u0;

    switch settings.initialEstimateMode
        case 'measured_angles'
            xhat0 = [yDev(1, 1); yDev(1, 2); 0; 0];
        case 'zero'
            xhat0 = zeros(4, 1);
        otherwise
            error('Unknown initialEstimateMode: %s', settings.initialEstimateMode);
    end

    [xhat, innovation] = kalman_run_filter(observerModel.sysData, Q, R, P0, xhat0, uDev, yDev);

    warmIdx = min(size(xhat, 1), max(1, floor(settings.warmupSeconds / settings.TsData) + 1));
    stateError = xhat(warmIdx:end, :) - ref.xDev(warmIdx:end, :);
    innovationWarm = innovation(warmIdx:end, :);

    angleNorm = stateError(:, 1:2) ./ normalization.angleScale(:).';
    rateNorm = stateError(:, 3:4) ./ normalization.rateScale(:).';

    score(i) = mean([angleNorm(:).^2; rateNorm(:).^2], 'omitnan');
    theta1Rmse(i) = sqrt(mean(stateError(:, 1).^2, 'omitnan'));
    theta2Rmse(i) = sqrt(mean(stateError(:, 2).^2, 'omitnan'));
    rate1Rmse(i) = sqrt(mean(stateError(:, 3).^2, 'omitnan'));
    rate2Rmse(i) = sqrt(mean(stateError(:, 4).^2, 'omitnan'));
    innovationRmse1(i) = sqrt(mean(innovationWarm(:, 1).^2, 'omitnan'));
    innovationRmse2(i) = sqrt(mean(innovationWarm(:, 2).^2, 'omitnan'));
end

experimentTable = table(names, score, theta1Rmse, theta2Rmse, rate1Rmse, rate2Rmse, ...
    innovationRmse1, innovationRmse2, ...
    'VariableNames', {'Experiment', 'Score', 'Theta1RMSE', 'Theta2RMSE', ...
    'Theta1DotRMSE', 'Theta2DotRMSE', 'InnovationTheta1RMSE', 'InnovationTheta2RMSE'});

metrics = struct();
metrics.experimentTable = experimentTable;
metrics.meanScore = mean(score, 'omitnan');
metrics.meanAngleRmse = mean([theta1Rmse; theta2Rmse], 'omitnan');
metrics.meanRateRmse = mean([rate1Rmse; rate2Rmse], 'omitnan');
metrics.meanInnovationRmse = mean([innovationRmse1; innovationRmse2], 'omitnan');
end

