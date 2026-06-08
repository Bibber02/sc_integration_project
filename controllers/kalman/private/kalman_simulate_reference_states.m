function references = kalman_simulate_reference_states(referenceModel, experiments, observerModel, settings)
%KALMAN_SIMULATE_REFERENCE_STATES Simulate reference angles and derive rates.

[z, initialStateValues] = kalman_make_iddata(experiments, settings);
modelForData = referenceModel;
modelForData = setinit(modelForData, 'Value', {
    initialStateValues(1, :);
    initialStateValues(2, :);
    initialStateValues(3, :);
    initialStateValues(4, :)
});
modelForData = setinit(modelForData, 'Fixed', {
    true(1, numel(experiments));
    true(1, numel(experiments));
    true(1, numel(experiments));
    true(1, numel(experiments))
});

simulatedOutput = [];
try
    opt = simOptions;
    opt.InitialCondition = 'm';
    simulatedOutput = sim(modelForData, z, opt);
catch ME
    warning('Reference model simulation failed; measured angles will be used as fallback. %s', ME.message);
end

references = struct( ...
    'name', {}, ...
    'xAbs', {}, ...
    'xDev', {}, ...
    'yAbs', {}, ...
    'yDev', {}, ...
    'uDev', {});

for i = 1:numel(experiments)
    yAbs = outputForExperiment(simulatedOutput, i);
    if isempty(yAbs) || any(~isfinite(yAbs), 'all')
        yAbs = experiments(i).y;
    end

    n = min(size(yAbs, 1), size(experiments(i).y, 1));
    yAbs = yAbs(1:n, :);
    dt = experiments(i).Ts;

    xAbs = zeros(n, 4);
    xAbs(:, 1:2) = yAbs;
    xAbs(:, 3) = gradient(yAbs(:, 1), dt);
    xAbs(:, 4) = gradient(yAbs(:, 2), dt);

    references(i).name = experiments(i).name;
    references(i).xAbs = xAbs;
    references(i).xDev = xAbs - observerModel.x0(:).';
    references(i).yAbs = xAbs(:, 1:2);
    references(i).yDev = xAbs(:, 1:2) - observerModel.y0(:).';
    references(i).uDev = experiments(i).u(1:n) - observerModel.u0;
end
end

function yAbs = outputForExperiment(simulatedOutput, index)
yAbs = [];
if isempty(simulatedOutput)
    return;
end

data = simulatedOutput.OutputData;
if iscell(data)
    if index <= numel(data)
        yAbs = double(data{index});
    end
else
    if index == 1
        yAbs = double(data);
    end
end
end
