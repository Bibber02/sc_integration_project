clear;
clc;

% ------------------------------------------------------------
% User settings
% ------------------------------------------------------------
% This script loads the parameter result from the full-system viscous
% grey-box identification result .mat file and creates the parameter vector p
% used by linearize_nonlinear_plant_viscous.m.
%
% Expected viscous full-system parameters:
%   p_a, p_b1, p_g1, p_u, p_0, p_b2, p_g2
%
% The coupling parameter p_c is then computed from p_g2:
%   p_c = (l1/g)*p_g2
%
% Final parameter vector:
%   p = [p_a; p_b1; p_g1; p_u; p_0; p_b2; p_g2; p_c]

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end

% Change this path to the .mat file saved by your full-system viscous ID script.
matFile = fullfile(scriptFolder, ...
    '..', '..', 'system_identification', 'full_system', 'grey_box', 'viscous', ...
    'full_system_id_viscous_result.mat');

% Measured distance from motor pivot to passive joint.
l1 = 0.10;   % [m]

% Gravitational acceleration.
g = 9.81;    % [m/s^2]

% Keep this true if you do not want to use a fitted constant torque/input bias
% in the linearized plant. This matches the old load_parameters.m behaviour,
% where p0 was overwritten by zero after loading the table.
forceP0ToZero = true;
manualP0Value = 0;

% ------------------------------------------------------------
% Load .mat file
% ------------------------------------------------------------

if ~isfile(matFile)
    error('Could not find viscous full-system result file:\n%s\nChange matFile at the top of load_parameters_viscous.m.', matFile);
end

S = load(matFile);
fieldNames = fieldnames(S);

% ------------------------------------------------------------
% Extract viscous parameters
% ------------------------------------------------------------

requiredNames = {'p_a', 'p_b1', 'p_g1', 'p_u', 'p_0', 'p_b2', 'p_g2'};
requiredKeys = strings(size(requiredNames));
for k = 1:numel(requiredNames)
    requiredKeys(k) = lower(regexprep(requiredNames{k}, '[^a-zA-Z0-9]', ''));
end

values = NaN(numel(requiredNames), 1);
parameterSource = '';

% First try: find a table with parameter names and values.
for f = 1:numel(fieldNames)

    candidate = S.(fieldNames{f});

    if ~istable(candidate)
        continue;
    end

    varNames = string(candidate.Properties.VariableNames);
    varKeys = lower(regexprep(varNames, '[^a-zA-Z0-9]', ''));

    nameCol = find(varKeys == "parameter" | varKeys == "name" | varKeys == "parameters", 1);
    valueCol = find(varKeys == "value" | varKeys == "values", 1);

    if isempty(valueCol) && width(candidate) >= 2
        valueCol = 2;
    end

    if isempty(nameCol) || isempty(valueCol)
        continue;
    end

    nameData = candidate{:, nameCol};
    valueData = candidate{:, valueCol};

    nameData = string(nameData);
    nameKeys = lower(regexprep(nameData, '[^a-zA-Z0-9]', ''));

    valuesTry = NaN(numel(requiredNames), 1);

    for k = 1:numel(requiredNames)
        hit = find(nameKeys == requiredKeys(k), 1);

        if ~isempty(hit)
            v = valueData(hit);

            if iscell(v)
                v = v{1};
            end

            if isstring(v) || ischar(v)
                v = str2double(v);
            end

            valuesTry(k) = double(v);
        end
    end

    if all(isfinite(valuesTry))
        values = valuesTry;
        parameterSource = sprintf('table "%s"', fieldNames{f});
        break;
    end
end

% Second try: find a 7-row table and use its second column in expected order.
if any(isnan(values))
    for f = 1:numel(fieldNames)

        candidate = S.(fieldNames{f});

        if istable(candidate) && height(candidate) == 7 && width(candidate) >= 2
            valuesTry = candidate{:, 2};

            if iscell(valuesTry)
                valuesTry = cellfun(@double, valuesTry);
            elseif isstring(valuesTry) || ischar(valuesTry)
                valuesTry = str2double(valuesTry);
            end

            valuesTry = double(valuesTry(:));

            if numel(valuesTry) == 7 && all(isfinite(valuesTry))
                values = valuesTry;
                parameterSource = sprintf('7-row table "%s", second column', fieldNames{f});
                break;
            end
        end
    end
end

% Third try: read parameters directly from a saved idnlgrey model object.
if any(isnan(values))
    possibleModelNames = {'modelEst', 'model_est', 'modelLocked', 'model_locked', 'modelFinal', 'model_final'};

    for m = 1:numel(possibleModelNames)

        modelName = possibleModelNames{m};

        if ~isfield(S, modelName)
            continue;
        end

        modelCandidate = S.(modelName);

        try
            params = modelCandidate.Parameters;
        catch
            continue;
        end

        paramNames = strings(numel(params), 1);
        paramValues = NaN(numel(params), 1);

        for j = 1:numel(params)
            paramNames(j) = string(params(j).Name);
            paramValues(j) = double(params(j).Value);
        end

        paramKeys = lower(regexprep(paramNames, '[^a-zA-Z0-9]', ''));
        valuesTry = NaN(numel(requiredNames), 1);

        for k = 1:numel(requiredNames)
            hit = find(paramKeys == requiredKeys(k), 1);
            if ~isempty(hit)
                valuesTry(k) = paramValues(hit);
            end
        end

        if all(isfinite(valuesTry))
            values = valuesTry;
            parameterSource = sprintf('model object "%s"', modelName);
            break;
        end
    end
end

if any(isnan(values))
    disp('Fields found in the loaded .mat file:');
    disp(fieldNames);
    error(['Could not extract all required viscous parameters. Required names are:\n' ...
           'p_a, p_b1, p_g1, p_u, p_0, p_b2, p_g2.\n' ...
           'Check that matFile points to the full-system viscous identification result.']);
end

fprintf('Loaded viscous full-system parameters from %s.\n', parameterSource);

% ------------------------------------------------------------
% Assign parameters
% ------------------------------------------------------------

p_a  = values(1);
p_b1 = values(2);
p_g1 = values(3);
p_u  = values(4);
p_0  = values(5);
p_b2 = values(6);
p_g2 = values(7);

if forceP0ToZero
    fprintf('p_0 loaded as %.6g, but overwritten by manualP0Value = %.6g.\n', p_0, manualP0Value);
    p_0 = manualP0Value;
end

% ------------------------------------------------------------
% Compute coupling parameter
% ------------------------------------------------------------
% In the reduced model:
%   p_c = gamma / beta = (l1/g)*p_g2

p_c = (l1/g)*p_g2;

% ------------------------------------------------------------
% Final parameter vector for linearization
% ------------------------------------------------------------

p = [
    p_a;
    p_b1;
    p_g1;
    p_u;
    p_0;
    p_b2;
    p_g2;
    p_c
];

fprintf('Created viscous parameter vector p with size %dx%d.\n', size(p,1), size(p,2));

parameterTable = array2table(p, ...
    'VariableNames', {'Value'}, ...
    'RowNames', {'p_a','p_b1','p_g1','p_u','p_0','p_b2','p_g2','p_c'});

disp(parameterTable);
