function [values, source] = kalman_reference_initial_parameters(projectPaths, parameterNames, fallbackValues)
%KALMAN_REFERENCE_INITIAL_PARAMETERS Warm-start reference ID when possible.

values = fallbackValues(:);
source = 'manual defaults plus passive-link result';

candidateFile = projectPaths.fullSystemStribeckResultWithP0;
if ~isfile(candidateFile)
    return;
end

S = load(candidateFile);
if ~isfield(S, 'resultTable') || ~istable(S.resultTable)
    return;
end

T = S.resultTable;
if ~all(ismember({'Parameter', 'Value'}, T.Properties.VariableNames))
    return;
end

tableNames = string(T.Parameter);

for k = 1:numel(parameterNames)
    idx = find(tableNames == string(parameterNames{k}), 1);
    if isempty(idx)
        continue;
    end

    value = T.Value(idx);
    if iscell(value)
        value = value{1};
    end
    value = double(value);

    if isfinite(value)
        values(k) = value;
    end
end

source = candidateFile;
end
