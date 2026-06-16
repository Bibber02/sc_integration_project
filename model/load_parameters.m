function [p, info] = load_parameters(settings)
%LOAD_PARAMETERS Load the selected no-p0 down-line identified model.
%
% Parameter vector convention used by nonlinearPlant, linearization and EKF:
%
%   p = [p_a; p_b1; p_c1; p_g1; p_u; ...
%        p_b2; p_g2; p_c2; p_sdelta2; v_s2; eps_v1; eps_v2; ...
%        theta_scale; theta1_offset; theta_abs_down_raw]
%
% There is deliberately no motor-bias p0 in this vector.
%
% Default selected result:
%   full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_sensor_calibration.mat
%   variable: modelStage3
%
% Change selectedStage below, or pass a settings struct, if you intentionally
% want another saved stage.

if nargin < 1 || isempty(settings)
    settings = struct();
elseif ischar(settings) || isstring(settings)
    settings = struct('selectedStage', char(settings));
end

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptFolder);

idFolder = fullfile(projectRoot, ...
    'system_identification', 'full_system', 'grey_box', 'stribeck');

selectedStage = settingOrDefault(settings, 'selectedStage', 'stage3_auto');

% Explicit override if needed:
%   load_parameters(struct('matFile', '...', 'modelVariable', 'modelStage4'))
if isfield(settings, 'matFile') && ~isempty(settings.matFile)
    matFile = char(settings.matFile);
    modelVariable = settingOrDefault(settings, 'modelVariable', '');
else
    [matFile, modelVariable] = resultFileForStage(idFolder, selectedStage);
end

% Fallback order. This keeps the control setup usable if you only have the
% LM files, or if you decide to use Stage 4 instead of Stage 3.
if ~isfile(matFile)
    candidates = {
        fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_sensor_calibration.mat'), 'modelStage3';
        fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage4_final_refinement.mat'), 'modelStage4';
        fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage3_sensor_calibration.mat'),      'modelStage3';
        fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage4_final_refinement.mat'),         'modelStage4';
        fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_stage3_sensor_calibration.mat'),                       'modelStage3';
        fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_stage4_final_refinement.mat'),                         'modelStage4';
        fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_80_20_stage3_sensor_calibration.mat'),                 'modelStage3';
        fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_80_20_stage4_final_refinement.mat'),                   'modelStage4';
        };

    found = false;
    for k = 1:size(candidates, 1)
        if isfile(candidates{k, 1})
            matFile = candidates{k, 1};
            modelVariable = candidates{k, 2};
            found = true;
            break;
        end
    end

    if ~found
        error(['Could not find a no-p0 down-line result .mat file in:\n  %s\n\n' ...
               'Expected for example:\n  %s\n\n' ...
               'Run the selected identification script first, or pass:\n' ...
               '  load_parameters(struct(''matFile'', ''your_file.mat'', ''modelVariable'', ''modelStage3''))'], ...
               idFolder, fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_sensor_calibration.mat'));
    end
end

S = load(matFile);

if isempty(modelVariable)
    modelVariable = findModelVariable(S);
end

parameterNames = {'p_a', 'p_b1', 'p_c1', 'p_g1', 'p_u', ...
    'p_b2', 'p_g2', 'p_c2', 'p_sdelta2', 'v_s2', 'eps_v1', 'eps_v2', ...
    'theta_scale', 'theta1_offset', 'theta_abs_down_raw'};

values = NaN(numel(parameterNames), 1);

if isfield(S, modelVariable)
    model = S.(modelVariable);
    for k = 1:numel(parameterNames)
        values(k) = getModelParameterValue(model, parameterNames{k});
    end
else
    % Fallback for result tables if the model object is not available.
    tableName = findResultTableVariable(S);
    if isempty(tableName)
        error('Could not find model variable "%s" or a result table in %s.', modelVariable, matFile);
    end
    T = S.(tableName);
    for k = 1:numel(parameterNames)
        values(k) = getTableParameterValue(T, parameterNames{k});
    end
end

if any(isnan(values))
    missing = parameterNames(isnan(values));
    error('Missing required no-p0 down-line parameters: %s', strjoin(missing, ', '));
end

p = values(:);

theta_scale = p(13);
theta1_offset = p(14);
theta_abs_down_raw = p(15);
theta2_offset = pi - theta_scale*theta_abs_down_raw - theta1_offset;

info = struct();
info.matFile = matFile;
info.modelVariable = modelVariable;
info.parameterNames = parameterNames;
info.theta2_offset = theta2_offset;
info.parameterVectorConvention = ['[p_a p_b1 p_c1 p_g1 p_u p_b2 p_g2 p_c2 ' ...
    'p_sdelta2 v_s2 eps_v1 eps_v2 theta_scale theta1_offset theta_abs_down_raw]'];

fprintf('Loaded no-p0 down-line parameter vector from:\n  %s\n', matFile);
fprintf('Using model variable: %s\n', modelVariable);
fprintf('Parameter vector length: %d\n', numel(p));
fprintf('theta_scale        = %.9g\n', theta_scale);
fprintf('theta1_offset      = %.9g rad = %.6g deg\n', theta1_offset, rad2deg(theta1_offset));
fprintf('theta_abs_down_raw = %.9g rad = %.6g deg\n', theta_abs_down_raw, rad2deg(theta_abs_down_raw));
fprintf('theta2_offset      = %.9g rad = %.6g deg   computed from down-line constraint\n', ...
    theta2_offset, rad2deg(theta2_offset));

end

%% Local helpers

function [matFile, modelVariable] = resultFileForStage(idFolder, selectedStage)
switch lower(char(selectedStage))
    case {'stage3_auto', 'auto_stage3', 'stage3'}
        matFile = fullfile(idFolder, ...
            'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_sensor_calibration.mat');
        modelVariable = 'modelStage3';

    case {'stage4_auto', 'auto_stage4', 'stage4'}
        matFile = fullfile(idFolder, ...
            'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage4_final_refinement.mat');
        modelVariable = 'modelStage4';

    case {'stage3_lm', 'lm_stage3'}
        matFile = fullfile(idFolder, ...
            'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage3_sensor_calibration.mat');
        modelVariable = 'modelStage3';

    case {'stage4_lm', 'lm_stage4'}
        matFile = fullfile(idFolder, ...
            'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage4_final_refinement.mat');
        modelVariable = 'modelStage4';

    otherwise
        error('Unknown selectedStage "%s". Use stage3_auto, stage4_auto, stage3_lm, or stage4_lm.', selectedStage);
end
end

function value = getModelParameterValue(model, parameterName)
value = NaN;
for k = 1:numel(model.Parameters)
    thisName = char(model.Parameters(k).Name);
    if strcmp(thisName, parameterName)
        rawValue = model.Parameters(k).Value;
        if iscell(rawValue)
            rawValue = rawValue{1};
        end
        value = double(rawValue);
        return;
    end
end
end

function value = getTableParameterValue(T, parameterName)
value = NaN;

if istable(T)
    % Most scripts store {'Parameter','Value',...}.
    if any(strcmp(T.Properties.VariableNames, 'Parameter')) && any(strcmp(T.Properties.VariableNames, 'Value'))
        names = T.Parameter;
        if iscell(names)
            mask = strcmp(names, parameterName);
        else
            mask = strcmp(cellstr(names), parameterName);
        end
        if any(mask)
            value = double(T.Value(find(mask, 1, 'first')));
            return;
        end
    end

    % Older generic resultTable had parameter names in column 1 and values
    % in column 2, or only values in column 2 in fixed known order.
    if width(T) >= 2
        try
            names = T{:, 1};
            if iscell(names) || isstring(names) || ischar(names)
                if ischar(names)
                    names = cellstr(names);
                end
                mask = strcmp(cellstr(names), parameterName);
                if any(mask)
                    value = double(T{find(mask, 1, 'first'), 2});
                    return;
                end
            end
        catch
        end
    end
end
end

function modelVariable = findModelVariable(S)
modelVariable = '';
preferred = {'modelStage3', 'modelStage4', 'modelStage2', 'modelEst'};
for k = 1:numel(preferred)
    if isfield(S, preferred{k})
        modelVariable = preferred{k};
        return;
    end
end
end

function tableName = findResultTableVariable(S)
tableName = '';
preferred = {'resultTableStage3', 'resultTableStage4', 'resultTableStage2', 'resultTable'};
for k = 1:numel(preferred)
    if isfield(S, preferred{k})
        tableName = preferred{k};
        return;
    end
end
end

function value = settingOrDefault(settings, fieldName, defaultValue)
if isfield(settings, fieldName) && ~isempty(settings.(fieldName))
    value = settings.(fieldName);
else
    value = defaultValue;
end
end
