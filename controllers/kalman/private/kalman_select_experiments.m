function selected = kalman_select_experiments(experiments, runNumbers)
%KALMAN_SELECT_EXPERIMENTS Select experiments by run number.

if isempty(experiments)
    selected = experiments;
    return;
end

mask = ismember([experiments.runNumber], runNumbers);
selected = experiments(mask);
end

