function settings = kalman_merge_settings(settings, overrides)
%KALMAN_MERGE_SETTINGS Shallow-merge user overrides into defaults.

fields = fieldnames(overrides);
for k = 1:numel(fields)
    settings.(fields{k}) = overrides.(fields{k});
end

if ~isfield(settings, 'referenceResultFile') || isempty(settings.referenceResultFile)
    settings.referenceResultFile = fullfile(settings.resultsFolder, 'kalman_reference_model.mat');
end
end

