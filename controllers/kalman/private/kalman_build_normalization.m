function normalization = kalman_build_normalization(references, R, settings)
%KALMAN_BUILD_NORMALIZATION Build scoring scales for angles and rates.

sensorStd = sqrt(max(diag(R), eps));
rateSamples = [];

for i = 1:numel(references)
    ref = references(i);
    warmIdx = findWarmupIndex(ref, settings);
    rateSamples = [rateSamples; ref.xDev(warmIdx:end, 3:4)]; %#ok<AGROW>
end

rateScale = sqrt(mean(rateSamples.^2, 1));
rateScale = max(rateScale(:), settings.minRateScale);

normalization = struct();
normalization.angleScale = sensorStd(:);
normalization.rateScale = rateScale(:);
end

function warmIdx = findWarmupIndex(ref, settings)
n = size(ref.xDev, 1);
warmIdx = min(n, max(1, floor(settings.warmupSeconds / settings.TsData) + 1));
end

