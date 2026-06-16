function [p, info] = load_parameters(settings)
%LOAD_PARAMETERS Load the selected no-p0 down-line identified model.
%
% This version deliberately avoids loading the complete identification .mat
% workspace. Those workspaces contain iddata/idnlgrey/idoptions objects, which
% can generate many warnings when they were saved with a different MATLAB /
% System Identification Toolbox version. For controller/EKF setup we only need
% the numeric parameter values, so this loader prefers the saved parameter CSV
% files and otherwise loads only the resultTable variable from the .mat file.
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
%   stage3_auto
%
% Examples:
%   p = load_parameters();
%   p = load_parameters('stage4_auto');
%   p = load_parameters('stage3_lm');
%   p = load_parameters(struct('parameterCsvFile','C:\path\stage3_parameters.csv'));

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

if isfield(settings, 'parameterCsvFile') && ~isempty(settings.parameterCsvFile)
    parameterCsvFile = char(settings.parameterCsvFile);
    matFile = '';
    tableVariable = '';
    modelVariable = '';
else
    if isfield(settings, 'matFile') && ~isempty(settings.matFile)
        matFile = char(settings.matFile);
        parameterCsvFile = settingOrDefault(settings, 'parameterCsvFile', '');
        modelVariable = settingOrDefault(settings, 'modelVariable', '');
        tableVariable = settingOrDefault(settings, 'tableVariable', '');
    else
        [matFile, parameterCsvFile, tableVariable, modelVariable] = resultFilesForStage(idFolder, selectedStage);
    end
end

% If the selected CSV/MAT is missing, fall back through the other reasonable
% no-p0 down-line stages. The CSV is preferred because it contains only plain
% numeric/table data and avoids toolbox object deserialization warnings.
if (isempty(parameterCsvFile) || ~isfile(parameterCsvFile)) && (isempty(matFile) || ~isfile(matFile))
    candidates = candidateFiles(idFolder);
    found = false;
    for k = 1:size(candidates, 1)
        csvTry = candidates{k, 2};
        matTry = candidates{k, 1};
        if isfile(csvTry) || isfile(matTry)
            matFile = matTry;
            parameterCsvFile = csvTry;
            tableVariable = candidates{k, 3};
            modelVariable = candidates{k, 4};
            found = true;
            break;
        end
    end

    if ~found
        error(['Could not find a no-p0 down-line parameter CSV or result .mat file in:\n  %s\n\n' ...
               'Expected for example:\n  %s\n\n' ...
               'Run the selected identification script first, or pass:\n' ...
               '  load_parameters(struct(''parameterCsvFile'', ''your_stage_parameters.csv''))'], ...
               idFolder, fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_parameters.csv'));
    end
end

parameterNames = {'p_a', 'p_b1', 'p_c1', 'p_g1', 'p_u', ...
    'p_b2', 'p_g2', 'p_c2', 'p_sdelta2', 'v_s2', 'eps_v1', 'eps_v2', ...
    'theta_scale', 'theta1_offset', 'theta_abs_down_raw'};

values = NaN(numel(parameterNames), 1);
sourceDescription = '';

% 1) Preferred path: read the saved parameter CSV.
if ~isempty(parameterCsvFile) && isfile(parameterCsvFile)
    T = readtable(parameterCsvFile);
    for k = 1:numel(parameterNames)
        values(k) = getTableParameterValue(T, parameterNames{k});
    end
    sourceDescription = parameterCsvFile;
end

% 2) Fallback: load only the result table variable from the MAT file.
% This avoids instantiating the many iddata/idnlgrey/idoptions objects saved
% in the workspace.
if any(isnan(values)) && ~isempty(matFile) && isfile(matFile)
    varsInFile = who('-file', matFile);

    tableCandidates = {};
    if ~isempty(tableVariable)
        tableCandidates{end+1} = tableVariable; %#ok<AGROW>
    end
    tableCandidates = [tableCandidates, {'resultTableStage3', 'resultTableStage4', 'resultTableStage2', 'resultTable'}]; %#ok<AGROW>

    loadedTable = false;
    for i = 1:numel(tableCandidates)
        thisTable = tableCandidates{i};
        if any(strcmp(varsInFile, thisTable))
            S = load(matFile, thisTable);
            T = S.(thisTable);
            for k = 1:numel(parameterNames)
                values(k) = getTableParameterValue(T, parameterNames{k});
            end
            sourceDescription = sprintf('%s : %s', matFile, thisTable);
            loadedTable = true;
            break;
        end
    end

    % 3) Last resort only: load the model object. This can reintroduce the
    % version warnings, but it is better than failing if no table exists.
    if any(isnan(values)) && ~loadedTable
        if isempty(modelVariable)
            modelVariable = findFirstExistingVariable(varsInFile, {'modelStage3', 'modelStage4', 'modelStage2', 'modelEst'});
        end
        if ~isempty(modelVariable) && any(strcmp(varsInFile, modelVariable))
            warning(['No parameter CSV/resultTable was found. Loading the idnlgrey model object instead. ' ...
                     'This may show MATLAB version warnings.']);
            S = load(matFile, modelVariable);
            model = S.(modelVariable);
            for k = 1:numel(parameterNames)
                values(k) = getModelParameterValue(model, parameterNames{k});
            end
            sourceDescription = sprintf('%s : %s', matFile, modelVariable);
        end
    end
end

if any(isnan(values))
    missing = parameterNames(isnan(values));
    error(['Missing required no-p0 down-line parameters: %s\n\n' ...
           'The selected result file was found, but the parameter table does not contain all required names.'], ...
           strjoin(missing, ', '));
end

p = values(:);

theta_scale = p(13);
theta1_offset = p(14);
theta_abs_down_raw = p(15);
theta2_offset = pi - theta_scale*theta_abs_down_raw - theta1_offset;

info = struct();
info.source = sourceDescription;
info.matFile = matFile;
info.parameterCsvFile = parameterCsvFile;
info.tableVariable = tableVariable;
info.modelVariable = modelVariable;
info.parameterNames = parameterNames;
info.theta2_offset = theta2_offset;
info.parameterVectorConvention = ['[p_a p_b1 p_c1 p_g1 p_u p_b2 p_g2 p_c2 ' ...
    'p_sdelta2 v_s2 eps_v1 eps_v2 theta_scale theta1_offset theta_abs_down_raw]'];

fprintf('Loaded no-p0 down-line parameter vector from:\n  %s\n', sourceDescription);
fprintf('Parameter vector length: %d\n', numel(p));
fprintf('theta_scale        = %.9g\n', theta_scale);
fprintf('theta1_offset      = %.9g rad = %.6g deg\n', theta1_offset, rad2deg(theta1_offset));
fprintf('theta_abs_down_raw = %.9g rad = %.6g deg\n', theta_abs_down_raw, rad2deg(theta_abs_down_raw));
fprintf('theta2_offset      = %.9g rad = %.6g deg   computed from down-line constraint\n', ...
    theta2_offset, rad2deg(theta2_offset));

end

%% Local helpers

function [matFile, csvFile, tableVariable, modelVariable] = resultFilesForStage(idFolder, selectedStage)
switch lower(char(selectedStage))
    case {'stage3_auto', 'auto_stage3', 'stage3'}
        base = 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3';
        matFile = fullfile(idFolder, [base '_sensor_calibration.mat']);
        csvFile = fullfile(idFolder, [base '_parameters.csv']);
        tableVariable = 'resultTableStage3';
        modelVariable = 'modelStage3';

    case {'stage4_auto', 'auto_stage4', 'stage4'}
        base = 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage4';
        matFile = fullfile(idFolder, [base '_final_refinement.mat']);
        csvFile = fullfile(idFolder, [base '_parameters.csv']);
        tableVariable = 'resultTableStage4';
        modelVariable = 'modelStage4';

    case {'stage2_auto', 'auto_stage2'}
        base = 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage2';
        matFile = fullfile(idFolder, [base '_active_offsets.mat']);
        csvFile = fullfile(idFolder, [base '_parameters.csv']);
        tableVariable = 'resultTableStage2';
        modelVariable = 'modelStage2';

    case {'stage3_lm', 'lm_stage3'}
        base = 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage3';
        matFile = fullfile(idFolder, [base '_sensor_calibration.mat']);
        csvFile = fullfile(idFolder, [base '_parameters.csv']);
        tableVariable = 'resultTableStage3';
        modelVariable = 'modelStage3';

    case {'stage4_lm', 'lm_stage4'}
        base = 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage4';
        matFile = fullfile(idFolder, [base '_final_refinement.mat']);
        csvFile = fullfile(idFolder, [base '_parameters.csv']);
        tableVariable = 'resultTableStage4';
        modelVariable = 'modelStage4';

    case {'stage2_lm', 'lm_stage2'}
        base = 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage2';
        matFile = fullfile(idFolder, [base '_active_offsets.mat']);
        csvFile = fullfile(idFolder, [base '_parameters.csv']);
        tableVariable = 'resultTableStage2';
        modelVariable = 'modelStage2';

    otherwise
        error('Unknown selectedStage "%s". Use stage3_auto, stage4_auto, stage3_lm, or stage4_lm.', selectedStage);
end
end

function candidates = candidateFiles(idFolder)
candidates = {
    fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_sensor_calibration.mat'), fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_parameters.csv'), 'resultTableStage3', 'modelStage3';
    fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage4_final_refinement.mat'),    fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage4_parameters.csv'), 'resultTableStage4', 'modelStage4';
    fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage3_sensor_calibration.mat'),      fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage3_parameters.csv'),      'resultTableStage3', 'modelStage3';
    fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage4_final_refinement.mat'),         fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_stage4_parameters.csv'),      'resultTableStage4', 'modelStage4';
    fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_80_20_stage3_sensor_calibration.mat'),                 fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_80_20_stage3_parameters.csv'),                 'resultTableStage3', 'modelStage3';
    fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_80_20_stage4_final_refinement.mat'),                   fullfile(idFolder, 'full_system_id_stribeck_calib_downline_no_p0_80_20_stage4_parameters.csv'),                 'resultTableStage4', 'modelStage4';
    };
end

function value = getTableParameterValue(T, parameterName)
value = NaN;

if ~istable(T)
    return;
end

if any(strcmp(T.Properties.VariableNames, 'Parameter')) && any(strcmp(T.Properties.VariableNames, 'Value'))
    names = T.Parameter;
    if iscell(names)
        mask = strcmp(names, parameterName);
    else
        mask = strcmp(cellstr(string(names)), parameterName);
    end
    if any(mask)
        rawValue = T.Value(find(mask, 1, 'first'));
        if iscell(rawValue)
            rawValue = rawValue{1};
        end
        value = double(rawValue);
        return;
    end
end

% Generic fallback: names in first column, values in second column.
if width(T) >= 2
    try
        names = T{:, 1};
        if ischar(names)
            names = cellstr(names);
        elseif isnumeric(names)
            return;
        else
            names = cellstr(string(names));
        end
        mask = strcmp(names, parameterName);
        if any(mask)
            rawValue = T{find(mask, 1, 'first'), 2};
            if iscell(rawValue)
                rawValue = rawValue{1};
            end
            value = double(rawValue);
            return;
        end
    catch
    end
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

function variableName = findFirstExistingVariable(varsInFile, preferred)
variableName = '';
for k = 1:numel(preferred)
    if any(strcmp(varsInFile, preferred{k}))
        variableName = preferred{k};
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
