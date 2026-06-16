clear;
clear functions;
clc;
close all;

%% ================================================================
% Full-system nonlinear grey-box identification, Stribeck friction version
%
% Changes relative to the earlier full-system script:
%   - the motor torque offset p_0 is removed from the parameter vector
%   - sensor calibration is represented by scale and offset corrections
%     on the already-calibrated measured angles:
%         theta_1_phys = theta_scale * theta_1_meas + theta1_offset
%         theta_2_phys = theta_scale * theta_2_meas + theta2_offset
%   - stage 1 estimates active parameters and calibration offsets with
%     p_c1 fixed and calibration scales fixed at 1
%   - stage 2 unlocks p_c1 and refines active parameters and offsets
%     with calibration scales still fixed at 1
%   - stage 2 results are saved before changing the calibration scales
%   - stage 3 fixes the dynamic parameters from stage 2 and estimates
%     the down-line constrained calibration correction parameters
%   - stage 4 performs a short final joint refinement of active dynamics
%     and down-line constrained calibration parameters with passive parameters fixed
%   - the experiment split is 80% training and 20% validation by run index
%
% Required model file on the MATLAB path:
%   greybox_id_full_stribeck_model_calib_downline_no_p0.m
% ================================================================

%% ================================================================
% Configuration
% ================================================================

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end

% This assumes the script is in:
%   full_system/grey_box/stribeck
% and the data folder is in:
%   full_system/measurement_data
% Change this single path if your folder structure is different.
dataFolder   = fullfile(scriptFolder, '..', '..', 'measurement_data');
prbsFolder   = fullfile(dataFolder, 'prbs');
chirpFolder  = fullfile(dataFolder, 'chirp');
modelsFolder = scriptFolder;
outputFolder = scriptFolder;

addpath(modelsFolder);

modelFile = 'greybox_id_full_stribeck_model_calib_downline_no_p0';

Ts = 0.01;
TsModel = 0;         % continuous-time grey-box model
restDurationForDownLine = 4.0;   % input starts at t = 4 s; use this section to define hanging-down raw line
inputSign = -1;      % original script used u = -u. Try +1 if the fit is mirrored.

amplitudes = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];

% 60/40 split by run index.
% Each selected run contributes both a PRBS and a chirp experiment.
% With 10 run indices and both input types enabled, this gives:
%   training:   6 run indices x 2 signal types = 12 experiments
%   validation: 4 run indices x 2 signal types =  8 experiments
%
% The validation runs are chosen as even-numbered runs spread over the amplitude range.
% The training set then keeps the remaining odd runs plus the highest amplitude.
idxVal = [2 4 6 8];
idxEst = setdiff(1:length(amplitudes), idxVal);

runPrbs = true;
runChirp = true;

% Hardware calibration values currently used in the real-time setup.
% These are not estimated directly because the saved data already contains
% calibrated angles, not raw ADC voltages. They are used only to print the
% equivalent corrected hardware gains and offsets after the fit.
a_theta_1_hw = 1.200374077783257;
b_theta_1_hw = 1.070535559513546;
a_theta_2_hw = 1.222725978474668;
b_theta_2_hw = 1.190099998448171;

% Initial guesses for calibration correction parameters applied to the
% already-calibrated measured angles.
theta_scale_init  = 1.0;
theta1_offset_init = 0.0;
theta_abs_down_raw_init = pi;  % overwritten after loading data

% Tight calibration bounds. Widen scale bounds only if the fitted values hit
% these limits and the result is still physically plausible.
thetaScaleMin = 0.97;
thetaScaleMax = 1.03;
thetaOffsetMin = -0.10;
thetaOffsetMax =  0.10;

useDiary = true;
consoleLogFile = 'full_system_id_stribeck_four_stage_calibration_downline_no_p0_60_40_theta2w15_auto_console_output.txt';

showInitialCompare = true;
showStage1Compare = true;
showStage2Compare = true;
showStage3Compare = true;
showStage4Compare = true;

figureUnits = 'pixels';
figurePosition = [80 80 1150 720];

% Passive-link parameters.
% The full-system fit should use the passive-link result directly, instead
% of copying numbers by hand. This keeps the full-system script consistent
% with the passive-link identification workflow.
passiveResultFolder = 'C:\School\Integration Project Systems and Control\sc_integration_project\system_identification\passive_link\grey_box\passive_link_id_3step_stribeck';
passiveResultFileName = 'passive_link_3step_stribeck_result_v5_no_helpers.mat';
passiveResultMatFile = fullfile(passiveResultFolder, passiveResultFileName);

% Fallback for copied/renamed result files.
if ~isfile(passiveResultMatFile)
    passiveResultFileName = 'passive_link_3step_stribeck_result.mat';
    passiveResultMatFile = fullfile(passiveResultFolder, passiveResultFileName);
end

if ~isfile(passiveResultMatFile)
    error('Passive-link result file not found in %s. Tried the v5 filename and passive_link_3step_stribeck_result.mat.', passiveResultFolder);
end

passiveLoaded = load(passiveResultMatFile);

if isfield(passiveLoaded, 'passiveID') && isfield(passiveLoaded.passiveID, 'mStrBest')
    passiveModel = passiveLoaded.passiveID.mStrBest;
elseif isfield(passiveLoaded, 'mStrBest')
    passiveModel = passiveLoaded.mStrBest;
else
    error('Passive-link result file does not contain passiveID.mStrBest or mStrBest.');
end

p_b2_val      = NaN;
p_g2_val      = NaN;
p_c2_val      = NaN;
p_sdelta2_val = NaN;
v_s2_val      = NaN;
eps_v2_val    = NaN;

for kPassivePar = 1:length(passiveModel.Parameters)
    passiveParName = char(passiveModel.Parameters(kPassivePar).Name);
    passiveParValue = passiveModel.Parameters(kPassivePar).Value;
    if iscell(passiveParValue)
        passiveParValue = passiveParValue{1};
    end
    passiveParValue = double(passiveParValue);

    switch passiveParName
        case 'p_b2'
            p_b2_val = passiveParValue;
        case 'p_g2'
            p_g2_val = passiveParValue;
        case 'p_c2'
            p_c2_val = passiveParValue;
        case {'p_sdelta2', 'p_sdelta'}
            p_sdelta2_val = passiveParValue;
        case {'v_s2', 'v_s'}
            v_s2_val = passiveParValue;
        case {'eps_v2', 'eps_v'}
            eps_v2_val = passiveParValue;
    end
end

if any(isnan([p_b2_val, p_g2_val, p_c2_val, p_sdelta2_val, v_s2_val, eps_v2_val]))
    error('Passive-link result file does not contain all required parameters.');
end

% Initial active/full-system parameter guesses.
p_a_init   = 14;
p_b1_init  = 670;
p_c1_init  = 0;
p_g1_init  = 97;
p_u_init   = 4000;
eps_v1_val = 0.05;

% Parameters:
%  1 p_a        active/full-system inertia ratio
%  2 p_b1       active joint effective viscous damping
%  3 p_c1       active joint Coulomb friction
%  4 p_g1       active-link gravity parameter
%  5 p_u        motor input gain
%  6 p_b2       passive joint viscous damping
%  7 p_g2       passive-link gravity parameter
%  8 p_c2       passive joint Coulomb friction
%  9 p_sdelta2  passive Stribeck static-minus-Coulomb term
% 10 v_s2       passive Stribeck velocity
% 11 eps_v1          active smoothing velocity
% 12 eps_v2          passive smoothing velocity
% 13 theta_scale        shared calibration scale correction
% 14 theta1_offset      theta_1 calibration offset correction
% 15 theta_abs_down_raw fixed raw value of theta_1 + theta_2 during initial rest
%    theta2_offset is computed from the measured initial-rest down-line constraint

parameterNames = {'p_a', 'p_b1', 'p_c1', 'p_g1', 'p_u', ...
    'p_b2', 'p_g2', 'p_c2', 'p_sdelta2', 'v_s2', 'eps_v1', 'eps_v2', ...
    'theta_scale', 'theta1_offset', 'theta_abs_down_raw'};

parameters = {
    p_a_init;
    p_b1_init;
    p_c1_init;
    p_g1_init;
    p_u_init;
    p_b2_val;
    p_g2_val;
    p_c2_val;
    p_sdelta2_val;
    v_s2_val;
    eps_v1_val;
    eps_v2_val;
    theta_scale_init;
    theta1_offset_init;
    theta_abs_down_raw_init
};

% Bounds for the Stribeck full-system model parameters.
% p_a is not tightly known from first-principles in this reduced model.
% Keep this bound configurable and inspect whether the fitted value still
% runs into the limit. If it does, the model is likely using p_a as a
% compensator for other modelling errors.
p_a_min = 5;
p_a_max = 150;

minimumValues = [p_a_min, 100, 0, 5, 500, 0, 20, 0, 0, 0.02, 0.001, 0.001, thetaScaleMin, thetaOffsetMin, pi-0.20];
maximumValues = [p_a_max, 4000, 800, 500, 30000, 20, 200, 20, 20, 20, 1, 1, thetaScaleMax, thetaOffsetMax, pi+0.20];

% Stage 1:
% Estimate p_a, p_b1, p_g1 and p_u.
% Keep p_c1 fixed at zero first, and keep all passive-link parameters fixed.
fixedStage1 = [false, false, true, false, false, ...
               true,  true,  true, true,  true, true, true, ...
               true,  false, true];

% Stage 2:
% Unlock p_c1 and estimate the final active/full-system parameter set.
% Keep passive-link parameters fixed.
fixedStage2 = [false, false, false, false, false, ...
               true,  true,  true,  true,  true, true, true, ...
               true,  false, true];

% Stage 3:
% Fix the dynamic parameters found in stage 2.
% Estimate the full calibration correction: scale and offset for both sensors.
% Passive parameters remain fixed because the previous passive-refit stage did
% not materially improve the validation fit.
fixedStage3 = [true, true, true, true, true, ...
               true, true, true, true, true, true, true, ...
               false, false, true];

% Stage 4:
% Short final joint refinement. Active/full-system parameters and calibration
% parameters are estimated together, while passive parameters and smoothing
% constants remain fixed. This lets the active dynamics adjust after the
% calibration scales have moved in stage 3, without reopening the passive
% parameter coupling.
fixedStage4 = [false, false, false, false, false, ...
               true,  true,  true,  true,  true, true, true, ...
               false, false, true];

stage2_p_c1_initial_value = 10;

% Estimation options.
% SearchMethod is set to 'auto'. With Optimization Toolbox available, MATLAB typically uses
% a bounded nonlinear least-squares solver for this grey-box prediction-error problem.
% The tolerances are deliberately not set extremely small. Once the step
% norm reaches around 1e-8, further improvement is not physically meaningful
% for this data/model combination. The try/catch assignments keep the script
% compatible with MATLAB releases that expose these tolerances under slightly
% different option names.
searchStepTolerance = 1e-8;
searchFunctionTolerance = 1e-6;

% Emphasise theta_2 more strongly in the prediction-error cost.
% Earlier scripts used diag([1, 8]); this version uses diag([1, 15]).
% Increase/decrease this single value if the theta_2 fit becomes over/under-prioritised.
theta2OutputWeight = 15;
outputWeightMatrix = diag([1, theta2OutputWeight]);

optStage1 = nlgreyestOptions;
optStage1.Display = 'Full';
optStage1.EstimateCovariance = true;
optStage1.SearchMethod = 'auto';
optStage1.SearchOptions.MaxIterations = 25;
optStage1.OutputWeight = outputWeightMatrix;
try, optStage1.SearchOptions.StepTolerance = searchStepTolerance; catch, end
try, optStage1.SearchOptions.FunctionTolerance = searchFunctionTolerance; catch, end
try, optStage1.SearchOptions.Advanced.TolX = searchStepTolerance; catch, end
try, optStage1.SearchOptions.Advanced.TolFun = searchFunctionTolerance; catch, end

optStage2 = nlgreyestOptions;
optStage2.Display = 'Full';
optStage2.EstimateCovariance = true;
optStage2.SearchMethod = 'auto';
optStage2.SearchOptions.MaxIterations = 35;
optStage2.OutputWeight = outputWeightMatrix;
try, optStage2.SearchOptions.StepTolerance = searchStepTolerance; catch, end
try, optStage2.SearchOptions.FunctionTolerance = searchFunctionTolerance; catch, end
try, optStage2.SearchOptions.Advanced.TolX = searchStepTolerance; catch, end
try, optStage2.SearchOptions.Advanced.TolFun = searchFunctionTolerance; catch, end

optStage3 = nlgreyestOptions;
optStage3.Display = 'Full';
optStage3.EstimateCovariance = true;
optStage3.SearchMethod = 'auto';
optStage3.SearchOptions.MaxIterations = 25;
optStage3.OutputWeight = outputWeightMatrix;
try, optStage3.SearchOptions.StepTolerance = searchStepTolerance; catch, end
try, optStage3.SearchOptions.FunctionTolerance = searchFunctionTolerance; catch, end
try, optStage3.SearchOptions.Advanced.TolX = searchStepTolerance; catch, end
try, optStage3.SearchOptions.Advanced.TolFun = searchFunctionTolerance; catch, end

optStage4 = nlgreyestOptions;
optStage4.Display = 'Full';
optStage4.EstimateCovariance = true;
optStage4.SearchMethod = 'auto';
optStage4.SearchOptions.MaxIterations = 20;
optStage4.OutputWeight = outputWeightMatrix;
try, optStage4.SearchOptions.StepTolerance = searchStepTolerance; catch, end
try, optStage4.SearchOptions.FunctionTolerance = searchFunctionTolerance; catch, end
try, optStage4.SearchOptions.Advanced.TolX = searchStepTolerance; catch, end
try, optStage4.SearchOptions.Advanced.TolFun = searchFunctionTolerance; catch, end

% Use the stored initial states for comparisons. This avoids hiding offset or
% model errors by re-estimating the initial condition during plotting.
compareOpt = compareOptions;
compareOpt.InitialCondition = 'model';

stage2MatFile = fullfile(outputFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage2_active_offsets.mat');
stage2ParameterCsvFile = fullfile(outputFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage2_parameters.csv');
stage2CorrelationCsvFile = fullfile(outputFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage2_correlations.csv');

stage3MatFile = fullfile(outputFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_sensor_calibration.mat');
stage3ParameterCsvFile = fullfile(outputFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_parameters.csv');
stage3CorrelationCsvFile = fullfile(outputFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage3_correlations.csv');

stage4MatFile = fullfile(outputFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage4_final_refinement.mat');
stage4ParameterCsvFile = fullfile(outputFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage4_parameters.csv');
stage4CorrelationCsvFile = fullfile(outputFolder, 'full_system_id_stribeck_calib_downline_no_p0_60_40_theta2w15_auto_stage4_correlations.csv');

%% ================================================================
% Start console logging
% ================================================================

if useDiary
    diary(fullfile(outputFolder, consoleLogFile));
end

fprintf('\n======================================================\n');
fprintf('Starting full-system Stribeck identification, no p_0\n');
fprintf('Sensor calibration scale/offset corrections are estimated as bounded grey-box parameters\n');
fprintf('theta2_offset is computed from theta_scale, theta1_offset, and the measured initial-rest down-line mean\n');
fprintf('======================================================\n\n');

fprintf('Data folder: %s\n', dataFolder);
fprintf('Model file:  %s\n', modelFile);
fprintf('theta_scale_init  = %.9g\n', theta_scale_init);
fprintf('theta1_offset_init = %.9g rad = %.6g deg\n', theta1_offset_init, rad2deg(theta1_offset_init));
fprintf('theta_abs_down_raw will be computed from the first %.2f s before input starts.\n\n', restDurationForDownLine);

fprintf('Output weighting: diag([1, %.6g])\n', theta2OutputWeight);
fprintf('Search method: auto\n');
fprintf('Identification runs: ');
fprintf('%d ', idxEst);
fprintf('\nValidation runs:     ');
fprintf('%d ', idxVal);
fprintf('\n\n');

modelPath = which(modelFile);
if isempty(modelPath)
    error(['MATLAB cannot find the model file "%s".\n' ...
           'Check modelsFolder and run: which %s'], modelFile, modelFile);
else
    fprintf('Using model file: %s\n\n', modelPath);
end

fprintf('Passive-link result file: %s\n', passiveResultMatFile);
fprintf('Initial passive-link parameters used in all stages:\n');
fprintf('  p_b2      = %.9g\n', p_b2_val);
fprintf('  p_g2      = %.9g\n', p_g2_val);
fprintf('  p_c2      = %.9g\n', p_c2_val);
fprintf('  p_sdelta2 = %.9g\n', p_sdelta2_val);
fprintf('  v_s2      = %.9g\n', v_s2_val);
fprintf('  eps_v2    = %.9g\n\n', eps_v2_val);

%% ================================================================
% Load data and build multi-experiment iddata objects
% ================================================================

zEst = [];
zVal = [];

x0Est = [];
x0Val = [];

experimentNamesEst = {};
experimentNamesVal = {};

nEst = 0;
nVal = 0;

thetaAbsDownRawEstSamples = [];
thetaAbsDownRawValSamples = [];

for kRun = 1:length(amplitudes)

    ampText = strrep(sprintf('%.2f', amplitudes(kRun)), '.', 'p');

    if runPrbs
        prbsFile = fullfile(prbsFolder, sprintf('fullsystem_prbs_A%s_run%02d.mat', ampText, kRun));

        if ~isfile(prbsFile)
            error('PRBS data file not found: %s', prbsFile);
        end

        load(prbsFile, 'theta_1', 'theta_2', 'u_ts');

        theta1 = double(squeeze(theta_1.Data(:)));
        theta2 = double(squeeze(theta_2.Data(:)));
        u = inputSign * double(squeeze(u_ts.Data(:)));

        nRestSamples = min(numel(theta1), max(1, round(restDurationForDownLine / Ts)));
        thetaAbsDownRawRest = theta1(1:nRestSamples) + theta2(1:nRestSamples);

        y = [theta1, theta2];

        z = iddata(y, u, Ts);
        z.Name = sprintf('prbs_run_%02d', kRun);
        z.InputName = {'u'};
        z.InputUnit = {'V'};
        z.OutputName = {'theta_1', 'theta_2'};
        z.OutputUnit = {'rad', 'rad'};
        z.TimeUnit = 's';

        x0 = [theta1(1); theta2(1); 0; 0];

        if ismember(kRun, idxEst)
            nEst = nEst + 1;

            if isempty(zEst)
                zEst = z;
            else
                zEst = merge(zEst, z);
            end

            x0Est(:, nEst) = x0;
            experimentNamesEst{nEst, 1} = z.Name;
            thetaAbsDownRawEstSamples = [thetaAbsDownRawEstSamples; thetaAbsDownRawRest(:)]; %#ok<AGROW>
        else
            nVal = nVal + 1;

            if isempty(zVal)
                zVal = z;
            else
                zVal = merge(zVal, z);
            end

            x0Val(:, nVal) = x0;
            experimentNamesVal{nVal, 1} = z.Name;
            thetaAbsDownRawValSamples = [thetaAbsDownRawValSamples; thetaAbsDownRawRest(:)]; %#ok<AGROW>
        end
    end

    if runChirp
        chirpFile = fullfile(chirpFolder, sprintf('fullsystem_chirp_A%s_run%02d.mat', ampText, kRun));

        if ~isfile(chirpFile)
            error('Chirp data file not found: %s', chirpFile);
        end

        load(chirpFile, 'theta_1', 'theta_2', 'u_ts');

        theta1 = double(squeeze(theta_1.Data(:)));
        theta2 = double(squeeze(theta_2.Data(:)));
        u = inputSign * double(squeeze(u_ts.Data(:)));

        nRestSamples = min(numel(theta1), max(1, round(restDurationForDownLine / Ts)));
        thetaAbsDownRawRest = theta1(1:nRestSamples) + theta2(1:nRestSamples);

        y = [theta1, theta2];

        z = iddata(y, u, Ts);
        z.Name = sprintf('chirp_run_%02d', kRun);
        z.InputName = {'u'};
        z.InputUnit = {'V'};
        z.OutputName = {'theta_1', 'theta_2'};
        z.OutputUnit = {'rad', 'rad'};
        z.TimeUnit = 's';

        x0 = [theta1(1); theta2(1); 0; 0];

        if ismember(kRun, idxEst)
            nEst = nEst + 1;

            if isempty(zEst)
                zEst = z;
            else
                zEst = merge(zEst, z);
            end

            x0Est(:, nEst) = x0;
            experimentNamesEst{nEst, 1} = z.Name;
            thetaAbsDownRawEstSamples = [thetaAbsDownRawEstSamples; thetaAbsDownRawRest(:)]; %#ok<AGROW>
        else
            nVal = nVal + 1;

            if isempty(zVal)
                zVal = z;
            else
                zVal = merge(zVal, z);
            end

            x0Val(:, nVal) = x0;
            experimentNamesVal{nVal, 1} = z.Name;
            thetaAbsDownRawValSamples = [thetaAbsDownRawValSamples; thetaAbsDownRawRest(:)]; %#ok<AGROW>
        end
    end
end

fprintf('Number of identification experiments: %d\n', nEst);
fprintf('Number of validation experiments:     %d\n\n', nVal);

if isempty(thetaAbsDownRawEstSamples)
    error('No identification rest samples were collected for the down-line constraint.');
end

theta_abs_down_raw = mean(thetaAbsDownRawEstSamples, 'omitnan');
theta_abs_down_raw_std = std(thetaAbsDownRawEstSamples, 'omitnan');
theta_abs_down_raw_val_mean = mean(thetaAbsDownRawValSamples, 'omitnan');
theta_abs_down_raw_val_std = std(thetaAbsDownRawValSamples, 'omitnan');

parameters{15} = theta_abs_down_raw;
minimumValues(15) = theta_abs_down_raw - 1e-9;
maximumValues(15) = theta_abs_down_raw + 1e-9;

fprintf('Down-line raw absolute angle from first %.2f s of training data:\n', restDurationForDownLine);
fprintf('  theta_abs_down_raw     = %.12g rad = %.8g deg\n', theta_abs_down_raw, rad2deg(theta_abs_down_raw));
fprintf('  training rest std      = %.6g rad = %.6g deg\n', theta_abs_down_raw_std, rad2deg(theta_abs_down_raw_std));
fprintf('  validation rest mean   = %.12g rad = %.8g deg\n', theta_abs_down_raw_val_mean, rad2deg(theta_abs_down_raw_val_mean));
fprintf('  validation rest std    = %.6g rad = %.6g deg\n\n', theta_abs_down_raw_val_std, rad2deg(theta_abs_down_raw_val_std));

if abs(theta_abs_down_raw - pi) > deg2rad(1)
    fprintf('NOTE: the measured initial-rest down-line is %.4g deg away from pi.\n', rad2deg(theta_abs_down_raw - pi));
    fprintf('      The constraint will map theta_1_raw + theta_2_raw = theta_abs_down_raw to physical pi.\n\n');
end

%% ================================================================
% Grey-box model setup
% ================================================================

order = [2 1 4];     % 2 outputs, 1 input, 4 states

initialStatesEst = {
    x0Est(1, :);
    x0Est(2, :);
    x0Est(3, :);
    x0Est(4, :)
};

model0 = idnlgrey(modelFile, order, parameters, initialStatesEst, TsModel);

for kPar = 1:length(parameterNames)
    model0.Parameters(kPar).Name = parameterNames{kPar};
    model0.Parameters(kPar).Minimum = minimumValues(kPar);
    model0.Parameters(kPar).Maximum = maximumValues(kPar);
    model0.Parameters(kPar).Fixed = fixedStage1(kPar);
end

model0.InputName = {'u'};
model0.InputUnit = {'V'};
model0.OutputName = {'theta_1', 'theta_2'};
model0.OutputUnit = {'rad', 'rad'};
model0.TimeUnit = 's';

model0.InitialStates(1).Name = 'theta_1_raw';
model0.InitialStates(1).Unit = 'rad';
model0.InitialStates(2).Name = 'theta_2_raw';
model0.InitialStates(2).Unit = 'rad';
model0.InitialStates(3).Name = 'theta_1_dot';
model0.InitialStates(3).Unit = 'rad/s';
model0.InitialStates(4).Name = 'theta_2_dot';
model0.InitialStates(4).Unit = 'rad/s';

% Keep the initial states fixed during nlgreyest, as in the original script.
model0 = setinit(model0, 'Fixed', {
    true(1, nEst);
    true(1, nEst);
    true(1, nEst);
    true(1, nEst)
});

%% ================================================================
% Compare initial model before estimation
% ================================================================

if showInitialCompare
    figure('Name', 'Initial model fit before estimation', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, model0, compareOpt);
    title('Initial model fit before estimation');
end

%% ================================================================
% Stage 1: estimate active parameters with p_c1 fixed
% ================================================================

fprintf('\n======================================================\n');
fprintf('Stage 1: estimating active parameters with p_c1 fixed\n');
fprintf('Free parameters: p_a, p_b1, p_g1, p_u, theta1_offset\n');
fprintf('Fixed calibration scales: theta_scale = 1\n');
fprintf('Fixed parameters: p_c1 and all passive-link parameters\n');
fprintf('======================================================\n\n');

modelStage1 = nlgreyest(zEst, model0, optStage1);
modelStage1.Name = 'Stage 1 active model, p_c1 fixed, calibration offsets estimated, no p_0';

fprintf('\nStage 1 parameter values:\n');
for kPar = 1:length(modelStage1.Parameters)
    if modelStage1.Parameters(kPar).Fixed
        fprintf('%-10s = %12.6g   fixed\n', modelStage1.Parameters(kPar).Name, modelStage1.Parameters(kPar).Value);
    else
        fprintf('%-10s = %12.6g   estimated\n', modelStage1.Parameters(kPar).Name, modelStage1.Parameters(kPar).Value);
    end
end

[~, fitStage1Est] = compare(zEst, modelStage1, compareOpt);

fitStage1Values = [];
if iscell(fitStage1Est)
    for kFit = 1:numel(fitStage1Est)
        fitStage1Values = [fitStage1Values; fitStage1Est{kFit}(:)]; %#ok<AGROW>
    end
else
    fitStage1Values = fitStage1Est(:);
end
fitStage1Values = fitStage1Values(isfinite(fitStage1Values));

fprintf('\nStage 1 mean identification fit: %.2f %%\n', mean(fitStage1Values));

if showStage1Compare
    figure('Name', 'Stage 1 identification fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, modelStage1, compareOpt);
    title('Stage 1 identification fit, p_{c1} fixed');
end

%% ================================================================
% Stage 2: unlock p_c1 and estimate final active parameter set
% ================================================================

fprintf('\n======================================================\n');
fprintf('Stage 2: estimating final active parameter set\n');
fprintf('Free parameters: p_a, p_b1, p_c1, p_g1, p_u, theta1_offset\n');
fprintf('Fixed calibration scales: theta_scale = 1\n');
fprintf('Fixed parameters: all passive-link parameters and smoothing constants\n');
fprintf('======================================================\n\n');

modelStage2Start = modelStage1;
modelStage2Start.Parameters(3).Value = stage2_p_c1_initial_value;

for kPar = 1:length(parameterNames)
    modelStage2Start.Parameters(kPar).Fixed = fixedStage2(kPar);
end

modelStage2 = nlgreyest(zEst, modelStage2Start, optStage2);
modelStage2.Name = 'Stage 2 full-system Stribeck model, calibration offsets estimated, no p_0, passive locked';

%% ================================================================
% Stage 2 validation and technical data
% ================================================================

initialStatesVal = {
    x0Val(1, :);
    x0Val(2, :);
    x0Val(3, :);
    x0Val(4, :)
};

parametersStage2Estimated = cell(length(parameterNames), 1);
for kPar = 1:length(parameterNames)
    parametersStage2Estimated{kPar} = modelStage2.Parameters(kPar).Value;
end

modelStage2Val = idnlgrey(modelFile, order, parametersStage2Estimated, initialStatesVal, TsModel);
for kPar = 1:length(parameterNames)
    modelStage2Val.Parameters(kPar).Name = parameterNames{kPar};
    modelStage2Val.Parameters(kPar).Minimum = minimumValues(kPar);
    modelStage2Val.Parameters(kPar).Maximum = maximumValues(kPar);
    modelStage2Val.Parameters(kPar).Fixed = true;
end
modelStage2Val.InputName = model0.InputName;
modelStage2Val.InputUnit = model0.InputUnit;
modelStage2Val.OutputName = model0.OutputName;
modelStage2Val.OutputUnit = model0.OutputUnit;
modelStage2Val.TimeUnit = model0.TimeUnit;
modelStage2Val.InitialStates(1).Name = model0.InitialStates(1).Name;
modelStage2Val.InitialStates(1).Unit = model0.InitialStates(1).Unit;
modelStage2Val.InitialStates(2).Name = model0.InitialStates(2).Name;
modelStage2Val.InitialStates(2).Unit = model0.InitialStates(2).Unit;
modelStage2Val.InitialStates(3).Name = model0.InitialStates(3).Name;
modelStage2Val.InitialStates(3).Unit = model0.InitialStates(3).Unit;
modelStage2Val.InitialStates(4).Name = model0.InitialStates(4).Name;
modelStage2Val.InitialStates(4).Unit = model0.InitialStates(4).Unit;
modelStage2Val = setinit(modelStage2Val, 'Fixed', {
    true(1, nVal);
    true(1, nVal);
    true(1, nVal);
    true(1, nVal)
});

[yStage2Est, fitStage2Est] = compare(zEst, modelStage2, compareOpt);
[yStage2Val, fitStage2Val] = compare(zVal, modelStage2Val, compareOpt);

fitStage2EstValues = [];
if iscell(fitStage2Est)
    for kFit = 1:numel(fitStage2Est)
        fitStage2EstValues = [fitStage2EstValues; fitStage2Est{kFit}(:)]; %#ok<AGROW>
    end
else
    fitStage2EstValues = fitStage2Est(:);
end
fitStage2EstValues = fitStage2EstValues(isfinite(fitStage2EstValues));

fitStage2ValValues = [];
if iscell(fitStage2Val)
    for kFit = 1:numel(fitStage2Val)
        fitStage2ValValues = [fitStage2ValValues; fitStage2Val{kFit}(:)]; %#ok<AGROW>
    end
else
    fitStage2ValValues = fitStage2Val(:);
end
fitStage2ValValues = fitStage2ValValues(isfinite(fitStage2ValValues));

fprintf('\nStage 2 identification fit, per experiment and output:\n');
disp(fitStage2Est);

fprintf('\nStage 2 validation fit, per experiment and output:\n');
disp(fitStage2Val);

fprintf('\nStage 2 mean identification fit: %.2f %%\n', mean(fitStage2EstValues));
fprintf('Stage 2 mean validation fit:     %.2f %%\n', mean(fitStage2ValValues));
theta2OffsetStage2 = pi - modelStage2.Parameters(13).Value * modelStage2.Parameters(15).Value - modelStage2.Parameters(14).Value;
fprintf('Stage 2 theta_scale:  %.9g   fixed\n', modelStage2.Parameters(13).Value);
fprintf('Stage 2 theta1_offset: %.9g rad = %.6g deg\n', modelStage2.Parameters(14).Value, rad2deg(modelStage2.Parameters(14).Value));
fprintf('Stage 2 theta_abs_down_raw:  %.9g   fixed\n', modelStage2.Parameters(15).Value);
fprintf('Stage 2 theta2_offset: %.9g rad = %.6g deg   computed from constraint\n', theta2OffsetStage2, rad2deg(theta2OffsetStage2));

if showStage2Compare
    figure('Name', 'Stage 2 identification fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, modelStage2, compareOpt);
    title('Stage 2 identification fit, active parameters estimated');

    figure('Name', 'Stage 2 validation fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zVal, modelStage2Val, compareOpt);
    title('Stage 2 validation fit, active parameters estimated');
end

fprintf('\n======================================================\n');
fprintf('Stage 2 parameter values and uncertainty\n');
fprintf('======================================================\n\n');

nParStage2 = length(modelStage2.Parameters);
parameterStage2 = cell(nParStage2, 1);
valueStage2 = zeros(nParStage2, 1);
fixedStage2Result = false(nParStage2, 1);
varianceStage2 = NaN(nParStage2, 1);
stdDevStage2 = NaN(nParStage2, 1);
relativeStdPercentStage2 = NaN(nParStage2, 1);

for kPar = 1:nParStage2
    parameterStage2{kPar} = modelStage2.Parameters(kPar).Name;
    valueStage2(kPar) = modelStage2.Parameters(kPar).Value;
    fixedStage2Result(kPar) = modelStage2.Parameters(kPar).Fixed;
end

freeIdxStage2 = find(~fixedStage2Result);
covFreeStage2 = getcov(modelStage2, 'value', 'free');

for kFree = 1:length(freeIdxStage2)
    varianceStage2(freeIdxStage2(kFree)) = covFreeStage2(kFree, kFree);
    stdDevStage2(freeIdxStage2(kFree)) = sqrt(max(covFreeStage2(kFree, kFree), 0));

    if abs(valueStage2(freeIdxStage2(kFree))) > eps
        relativeStdPercentStage2(freeIdxStage2(kFree)) = 100 * stdDevStage2(freeIdxStage2(kFree)) / abs(valueStage2(freeIdxStage2(kFree)));
    end
end

resultTableStage2 = table(parameterStage2, valueStage2, fixedStage2Result, varianceStage2, stdDevStage2, relativeStdPercentStage2, ...
    'VariableNames', {'Parameter', 'Value', 'Fixed', 'Variance', 'StdDev', 'RelativeStdPercent'});

disp(resultTableStage2);

fprintf('\nStage 2 free-parameter covariance matrix:\n');
disp(covFreeStage2);

freeNamesStage2 = parameterStage2(freeIdxStage2);
stdFreeStage2 = sqrt(max(diag(covFreeStage2), 0));
corrFreeStage2 = covFreeStage2 ./ (stdFreeStage2 * stdFreeStage2.');
corrFreeStage2(1:size(corrFreeStage2, 1)+1:end) = 1;

corrVarNamesStage2 = matlab.lang.makeValidName(freeNamesStage2);
corrTableStage2 = array2table(corrFreeStage2, 'VariableNames', corrVarNamesStage2, 'RowNames', freeNamesStage2);

fprintf('\nStage 2 free-parameter correlation matrix:\n');
disp(corrTableStage2);

fprintf('\nStage 2 loss function: %g\n', modelStage2.Report.Fit.LossFcn);

stage2SaveTime = datetime('now');
save(stage2MatFile);
writetable(resultTableStage2, stage2ParameterCsvFile);
writetable(corrTableStage2, stage2CorrelationCsvFile, 'WriteRowNames', true);

fprintf('\nSaved Stage 2 workspace and technical data before Stage 3:\n');
fprintf('  %s\n', stage2MatFile);
fprintf('  %s\n', stage2ParameterCsvFile);
fprintf('  %s\n\n', stage2CorrelationCsvFile);

%% ================================================================
% Stage 3: fix dynamics and estimate sensor calibration parameters
% ================================================================

fprintf('\n======================================================\n');
fprintf('Stage 3: estimating sensor calibration scale and offset corrections\n');
fprintf('Fixed parameters: active parameters, passive parameters, smoothing constants\n');
fprintf('Free parameters: theta_scale, theta1_offset\n');
fprintf('======================================================\n\n');

modelStage3Start = modelStage2;

for kPar = 1:length(parameterNames)
    modelStage3Start.Parameters(kPar).Fixed = fixedStage3(kPar);
end

% Stage 3 has only two free calibration parameters. A small deterministic
% multistart grid is cheap and reduces the risk that the calibration stage
% gets stuck in a poor local minimum after the active dynamics are fixed.
thetaScaleStarts = [0.99, 1.00, 1.01];
theta1OffsetStarts = modelStage2.Parameters(14).Value + [-0.03, 0, 0.03];
thetaScaleStarts = unique(max(thetaScaleMin, min(thetaScaleMax, thetaScaleStarts)));
theta1OffsetStarts = unique(max(thetaOffsetMin, min(thetaOffsetMax, theta1OffsetStarts)));

bestStage3Loss = Inf;
modelStage3 = [];
stage3StartResults = [];
startCounter = 0;

for kScaleStart = 1:numel(thetaScaleStarts)
    for kOffsetStart = 1:numel(theta1OffsetStarts)
        startCounter = startCounter + 1;
        modelStage3TryStart = modelStage3Start;
        modelStage3TryStart.Parameters(13).Value = thetaScaleStarts(kScaleStart);
        modelStage3TryStart.Parameters(14).Value = theta1OffsetStarts(kOffsetStart);

        fprintf('Stage 3 calibration start %d: theta_scale = %.6g, theta1_offset = %.6g rad\n', ...
            startCounter, modelStage3TryStart.Parameters(13).Value, modelStage3TryStart.Parameters(14).Value);

        try
            modelStage3Try = nlgreyest(zEst, modelStage3TryStart, optStage3);
            currentLoss = modelStage3Try.Report.Fit.LossFcn;
            stage3StartResults = [stage3StartResults; ... %#ok<AGROW>
                startCounter, modelStage3TryStart.Parameters(13).Value, modelStage3TryStart.Parameters(14).Value, ...
                modelStage3Try.Parameters(13).Value, modelStage3Try.Parameters(14).Value, currentLoss];

            fprintf('  completed: theta_scale = %.9g, theta1_offset = %.9g, loss = %.9g\n', ...
                modelStage3Try.Parameters(13).Value, modelStage3Try.Parameters(14).Value, currentLoss);

            if currentLoss < bestStage3Loss
                bestStage3Loss = currentLoss;
                modelStage3 = modelStage3Try;
            end
        catch ME
            warning('Stage 3 calibration start %d failed: %s', startCounter, ME.message);
            stage3StartResults = [stage3StartResults; ... %#ok<AGROW>
                startCounter, modelStage3TryStart.Parameters(13).Value, modelStage3TryStart.Parameters(14).Value, NaN, NaN, NaN];
        end
    end
end

if isempty(modelStage3)
    error('All Stage 3 calibration starts failed.');
end

stage3StartResultsTable = array2table(stage3StartResults, ...
    'VariableNames', {'Start', 'InitialThetaScale', 'InitialTheta1Offset', ...
    'FinalThetaScale', 'FinalTheta1Offset', 'LossFcn'});

fprintf('\nStage 3 calibration multistart results:\n');
disp(stage3StartResultsTable);

modelStage3.Name = 'Stage 3 full-system Stribeck model, dynamics fixed, down-line calibration refined';

%% ================================================================
% Stage 3 validation and technical data
% ================================================================

parametersStage3Estimated = cell(length(parameterNames), 1);
for kPar = 1:length(parameterNames)
    parametersStage3Estimated{kPar} = modelStage3.Parameters(kPar).Value;
end

modelStage3Val = idnlgrey(modelFile, order, parametersStage3Estimated, initialStatesVal, TsModel);
for kPar = 1:length(parameterNames)
    modelStage3Val.Parameters(kPar).Name = parameterNames{kPar};
    modelStage3Val.Parameters(kPar).Minimum = minimumValues(kPar);
    modelStage3Val.Parameters(kPar).Maximum = maximumValues(kPar);
    modelStage3Val.Parameters(kPar).Fixed = true;
end
modelStage3Val.InputName = model0.InputName;
modelStage3Val.InputUnit = model0.InputUnit;
modelStage3Val.OutputName = model0.OutputName;
modelStage3Val.OutputUnit = model0.OutputUnit;
modelStage3Val.TimeUnit = model0.TimeUnit;
modelStage3Val.InitialStates(1).Name = model0.InitialStates(1).Name;
modelStage3Val.InitialStates(1).Unit = model0.InitialStates(1).Unit;
modelStage3Val.InitialStates(2).Name = model0.InitialStates(2).Name;
modelStage3Val.InitialStates(2).Unit = model0.InitialStates(2).Unit;
modelStage3Val.InitialStates(3).Name = model0.InitialStates(3).Name;
modelStage3Val.InitialStates(3).Unit = model0.InitialStates(3).Unit;
modelStage3Val.InitialStates(4).Name = model0.InitialStates(4).Name;
modelStage3Val.InitialStates(4).Unit = model0.InitialStates(4).Unit;
modelStage3Val = setinit(modelStage3Val, 'Fixed', {
    true(1, nVal);
    true(1, nVal);
    true(1, nVal);
    true(1, nVal)
});

[yStage3Est, fitStage3Est] = compare(zEst, modelStage3, compareOpt);
[yStage3Val, fitStage3Val] = compare(zVal, modelStage3Val, compareOpt);

fitStage3EstValues = [];
if iscell(fitStage3Est)
    for kFit = 1:numel(fitStage3Est)
        fitStage3EstValues = [fitStage3EstValues; fitStage3Est{kFit}(:)]; %#ok<AGROW>
    end
else
    fitStage3EstValues = fitStage3Est(:);
end
fitStage3EstValues = fitStage3EstValues(isfinite(fitStage3EstValues));

fitStage3ValValues = [];
if iscell(fitStage3Val)
    for kFit = 1:numel(fitStage3Val)
        fitStage3ValValues = [fitStage3ValValues; fitStage3Val{kFit}(:)]; %#ok<AGROW>
    end
else
    fitStage3ValValues = fitStage3Val(:);
end
fitStage3ValValues = fitStage3ValValues(isfinite(fitStage3ValValues));

fprintf('\nStage 3 identification fit, per experiment and output:\n');
disp(fitStage3Est);

fprintf('\nStage 3 validation fit, per experiment and output:\n');
disp(fitStage3Val);

fprintf('\nStage 3 mean identification fit: %.2f %%\n', mean(fitStage3EstValues));
fprintf('Stage 3 mean validation fit:     %.2f %%\n', mean(fitStage3ValValues));
theta2OffsetStage3 = pi - modelStage3.Parameters(13).Value * modelStage3.Parameters(15).Value - modelStage3.Parameters(14).Value;
fprintf('Stage 3 theta_scale:  %.9g\n', modelStage3.Parameters(13).Value);
fprintf('Stage 3 theta1_offset: %.9g rad = %.6g deg\n', modelStage3.Parameters(14).Value, rad2deg(modelStage3.Parameters(14).Value));
fprintf('Stage 3 theta_abs_down_raw:  %.9g\n', modelStage3.Parameters(15).Value);
fprintf('Stage 3 theta2_offset: %.9g rad = %.6g deg   computed from constraint\n', theta2OffsetStage3, rad2deg(theta2OffsetStage3));

a_theta_1_corrected = modelStage3.Parameters(13).Value * a_theta_1_hw;
b_theta_1_corrected = b_theta_1_hw - modelStage3.Parameters(14).Value / a_theta_1_corrected;
a_theta_2_corrected = modelStage3.Parameters(13).Value * a_theta_2_hw;
b_theta_2_corrected = b_theta_2_hw - theta2OffsetStage3 / a_theta_2_corrected;

fprintf('Equivalent hardware calibration after Stage 3:\n');
fprintf('  a_theta_1 = %.15g\n', a_theta_1_corrected);
fprintf('  b_theta_1 = %.15g\n', b_theta_1_corrected);
fprintf('  a_theta_2 = %.15g\n', a_theta_2_corrected);
fprintf('  b_theta_2 = %.15g\n', b_theta_2_corrected);


if showStage3Compare
    figure('Name', 'Stage 3 identification fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, modelStage3, compareOpt);
    title('Stage 3 identification fit, sensor calibration refined');

    figure('Name', 'Stage 3 validation fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zVal, modelStage3Val, compareOpt);
    title('Stage 3 validation fit, sensor calibration refined');
end

fprintf('\n======================================================\n');
fprintf('Stage 3 parameter values and uncertainty\n');
fprintf('======================================================\n\n');

nParStage3 = length(modelStage3.Parameters);
parameterStage3 = cell(nParStage3, 1);
valueStage3 = zeros(nParStage3, 1);
fixedStage3Result = false(nParStage3, 1);
varianceStage3 = NaN(nParStage3, 1);
stdDevStage3 = NaN(nParStage3, 1);
relativeStdPercentStage3 = NaN(nParStage3, 1);

for kPar = 1:nParStage3
    parameterStage3{kPar} = modelStage3.Parameters(kPar).Name;
    valueStage3(kPar) = modelStage3.Parameters(kPar).Value;
    fixedStage3Result(kPar) = modelStage3.Parameters(kPar).Fixed;
end

freeIdxStage3 = find(~fixedStage3Result);
covFreeStage3 = getcov(modelStage3, 'value', 'free');

for kFree = 1:length(freeIdxStage3)
    varianceStage3(freeIdxStage3(kFree)) = covFreeStage3(kFree, kFree);
    stdDevStage3(freeIdxStage3(kFree)) = sqrt(max(covFreeStage3(kFree, kFree), 0));

    if abs(valueStage3(freeIdxStage3(kFree))) > eps
        relativeStdPercentStage3(freeIdxStage3(kFree)) = 100 * stdDevStage3(freeIdxStage3(kFree)) / abs(valueStage3(freeIdxStage3(kFree)));
    end
end

resultTableStage3 = table(parameterStage3, valueStage3, fixedStage3Result, varianceStage3, stdDevStage3, relativeStdPercentStage3, ...
    'VariableNames', {'Parameter', 'Value', 'Fixed', 'Variance', 'StdDev', 'RelativeStdPercent'});

disp(resultTableStage3);

fprintf('\nStage 3 free-parameter covariance matrix:\n');
disp(covFreeStage3);

freeNamesStage3 = parameterStage3(freeIdxStage3);
stdFreeStage3 = sqrt(max(diag(covFreeStage3), 0));
corrFreeStage3 = covFreeStage3 ./ (stdFreeStage3 * stdFreeStage3.');
corrFreeStage3(1:size(corrFreeStage3, 1)+1:end) = 1;

corrVarNamesStage3 = matlab.lang.makeValidName(freeNamesStage3);
corrTableStage3 = array2table(corrFreeStage3, 'VariableNames', corrVarNamesStage3, 'RowNames', freeNamesStage3);

fprintf('\nStage 3 free-parameter correlation matrix:\n');
disp(corrTableStage3);

fprintf('\nStage 3 loss function: %g\n', modelStage3.Report.Fit.LossFcn);

% Compatibility aliases for later scripts that expect the old names.
modelEst = modelStage3;
modelVal = modelStage3Val;
yEst = yStage3Est;
yVal = yStage3Val;
fitEst = fitStage3Est;
fitVal = fitStage3Val;
resultTable = resultTableStage3;
covFree = covFreeStage3;
corrTable = corrTableStage3;

stage3SaveTime = datetime('now');
save(stage3MatFile);
writetable(resultTableStage3, stage3ParameterCsvFile);
writetable(corrTableStage3, stage3CorrelationCsvFile, 'WriteRowNames', true);

fprintf('\nSaved Stage 3 workspace and technical data:\n');
fprintf('  %s\n', stage3MatFile);
fprintf('  %s\n', stage3ParameterCsvFile);

fprintf('  %s\n', stage3CorrelationCsvFile);

%% ================================================================
% Stage 4: short final joint refinement of active dynamics and calibration
% ================================================================

fprintf('\n======================================================\n');
fprintf('Stage 4: short final joint refinement, diagnostic only\n');
fprintf('Free parameters: active parameters and down-line constrained calibration corrections\n');
fprintf('Fixed parameters: passive parameters and smoothing constants\n');
fprintf('======================================================\n\n');

modelStage4Start = modelStage3;

for kPar = 1:length(parameterNames)
    modelStage4Start.Parameters(kPar).Fixed = fixedStage4(kPar);
end

modelStage4 = nlgreyest(zEst, modelStage4Start, optStage4);
modelStage4.Name = 'Stage 4 full-system Stribeck model, active dynamics and sensor calibration jointly refined';

%% ================================================================
% Stage 4 validation and technical data
% ================================================================

parametersStage4Estimated = cell(length(parameterNames), 1);
for kPar = 1:length(parameterNames)
    parametersStage4Estimated{kPar} = modelStage4.Parameters(kPar).Value;
end

modelStage4Val = idnlgrey(modelFile, order, parametersStage4Estimated, initialStatesVal, TsModel);
for kPar = 1:length(parameterNames)
    modelStage4Val.Parameters(kPar).Name = parameterNames{kPar};
    modelStage4Val.Parameters(kPar).Minimum = minimumValues(kPar);
    modelStage4Val.Parameters(kPar).Maximum = maximumValues(kPar);
    modelStage4Val.Parameters(kPar).Fixed = true;
end
modelStage4Val.InputName = model0.InputName;
modelStage4Val.InputUnit = model0.InputUnit;
modelStage4Val.OutputName = model0.OutputName;
modelStage4Val.OutputUnit = model0.OutputUnit;
modelStage4Val.TimeUnit = model0.TimeUnit;
modelStage4Val.InitialStates(1).Name = model0.InitialStates(1).Name;
modelStage4Val.InitialStates(1).Unit = model0.InitialStates(1).Unit;
modelStage4Val.InitialStates(2).Name = model0.InitialStates(2).Name;
modelStage4Val.InitialStates(2).Unit = model0.InitialStates(2).Unit;
modelStage4Val.InitialStates(3).Name = model0.InitialStates(3).Name;
modelStage4Val.InitialStates(3).Unit = model0.InitialStates(3).Unit;
modelStage4Val.InitialStates(4).Name = model0.InitialStates(4).Name;
modelStage4Val.InitialStates(4).Unit = model0.InitialStates(4).Unit;
modelStage4Val = setinit(modelStage4Val, 'Fixed', {
    true(1, nVal);
    true(1, nVal);
    true(1, nVal);
    true(1, nVal)
});

[yStage4Est, fitStage4Est] = compare(zEst, modelStage4, compareOpt);
[yStage4Val, fitStage4Val] = compare(zVal, modelStage4Val, compareOpt);

fitStage4EstValues = [];
if iscell(fitStage4Est)
    for kFit = 1:numel(fitStage4Est)
        fitStage4EstValues = [fitStage4EstValues; fitStage4Est{kFit}(:)]; %#ok<AGROW>
    end
else
    fitStage4EstValues = fitStage4Est(:);
end
fitStage4EstValues = fitStage4EstValues(isfinite(fitStage4EstValues));

fitStage4ValValues = [];
if iscell(fitStage4Val)
    for kFit = 1:numel(fitStage4Val)
        fitStage4ValValues = [fitStage4ValValues; fitStage4Val{kFit}(:)]; %#ok<AGROW>
    end
else
    fitStage4ValValues = fitStage4Val(:);
end
fitStage4ValValues = fitStage4ValValues(isfinite(fitStage4ValValues));

fprintf('\nStage 4 identification fit, per experiment and output:\n');
disp(fitStage4Est);

fprintf('\nStage 4 validation fit, per experiment and output:\n');
disp(fitStage4Val);

fprintf('\nStage 4 mean identification fit: %.2f %%\n', mean(fitStage4EstValues));
fprintf('Stage 4 mean validation fit:     %.2f %%\n', mean(fitStage4ValValues));
theta2OffsetStage4 = pi - modelStage4.Parameters(13).Value * modelStage4.Parameters(15).Value - modelStage4.Parameters(14).Value;
fprintf('Stage 4 theta_scale:  %.9g\n', modelStage4.Parameters(13).Value);
fprintf('Stage 4 theta1_offset: %.9g rad = %.6g deg\n', modelStage4.Parameters(14).Value, rad2deg(modelStage4.Parameters(14).Value));
fprintf('Stage 4 theta_abs_down_raw:  %.9g\n', modelStage4.Parameters(15).Value);
fprintf('Stage 4 theta2_offset: %.9g rad = %.6g deg   computed from constraint\n', theta2OffsetStage4, rad2deg(theta2OffsetStage4));

a_theta_1_corrected_stage4 = modelStage4.Parameters(13).Value * a_theta_1_hw;
b_theta_1_corrected_stage4 = b_theta_1_hw - modelStage4.Parameters(14).Value / a_theta_1_corrected_stage4;
a_theta_2_corrected_stage4 = modelStage4.Parameters(13).Value * a_theta_2_hw;
b_theta_2_corrected_stage4 = b_theta_2_hw - theta2OffsetStage4 / a_theta_2_corrected_stage4;

fprintf('Equivalent hardware calibration after Stage 4:\n');
fprintf('  a_theta_1 = %.15g\n', a_theta_1_corrected_stage4);
fprintf('  b_theta_1 = %.15g\n', b_theta_1_corrected_stage4);
fprintf('  a_theta_2 = %.15g\n', a_theta_2_corrected_stage4);
fprintf('  b_theta_2 = %.15g\n', b_theta_2_corrected_stage4);

if showStage4Compare
    figure('Name', 'Stage 4 identification fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zEst, modelStage4, compareOpt);
    title('Stage 4 identification fit, final joint refinement');

    figure('Name', 'Stage 4 validation fit', ...
        'Units', figureUnits, 'Position', figurePosition);
    compare(zVal, modelStage4Val, compareOpt);
    title('Stage 4 validation fit, final joint refinement');
end

fprintf('\n======================================================\n');
fprintf('Stage 4 parameter values and uncertainty\n');
fprintf('======================================================\n\n');

nParStage4 = length(modelStage4.Parameters);
parameterStage4 = cell(nParStage4, 1);
valueStage4 = zeros(nParStage4, 1);
fixedStage4Result = false(nParStage4, 1);
varianceStage4 = NaN(nParStage4, 1);
stdDevStage4 = NaN(nParStage4, 1);
relativeStdPercentStage4 = NaN(nParStage4, 1);

for kPar = 1:nParStage4
    parameterStage4{kPar} = modelStage4.Parameters(kPar).Name;
    valueStage4(kPar) = modelStage4.Parameters(kPar).Value;
    fixedStage4Result(kPar) = modelStage4.Parameters(kPar).Fixed;
end

freeIdxStage4 = find(~fixedStage4Result);
covFreeStage4 = getcov(modelStage4, 'value', 'free');

for kFree = 1:length(freeIdxStage4)
    varianceStage4(freeIdxStage4(kFree)) = covFreeStage4(kFree, kFree);
    stdDevStage4(freeIdxStage4(kFree)) = sqrt(max(covFreeStage4(kFree, kFree), 0));

    if abs(valueStage4(freeIdxStage4(kFree))) > eps
        relativeStdPercentStage4(freeIdxStage4(kFree)) = 100 * stdDevStage4(freeIdxStage4(kFree)) / abs(valueStage4(freeIdxStage4(kFree)));
    end
end

resultTableStage4 = table(parameterStage4, valueStage4, fixedStage4Result, varianceStage4, stdDevStage4, relativeStdPercentStage4, ...
    'VariableNames', {'Parameter', 'Value', 'Fixed', 'Variance', 'StdDev', 'RelativeStdPercent'});

disp(resultTableStage4);

fprintf('\nStage 4 free-parameter covariance matrix:\n');
disp(covFreeStage4);

freeNamesStage4 = parameterStage4(freeIdxStage4);
stdFreeStage4 = sqrt(max(diag(covFreeStage4), 0));
corrFreeStage4 = covFreeStage4 ./ (stdFreeStage4 * stdFreeStage4.');
corrFreeStage4(1:size(corrFreeStage4, 1)+1:end) = 1;

corrVarNamesStage4 = matlab.lang.makeValidName(freeNamesStage4);
corrTableStage4 = array2table(corrFreeStage4, 'VariableNames', corrVarNamesStage4, 'RowNames', freeNamesStage4);

fprintf('\nStage 4 free-parameter correlation matrix:\n');
disp(corrTableStage4);

fprintf('\nStage 4 loss function: %g\n', modelStage4.Report.Fit.LossFcn);

% Compatibility aliases for later scripts that expect the old final-model names.
modelEst = modelStage4;
modelVal = modelStage4Val;
yEst = yStage4Est;
yVal = yStage4Val;
fitEst = fitStage4Est;
fitVal = fitStage4Val;
resultTable = resultTableStage4;
covFree = covFreeStage4;
corrTable = corrTableStage4;

stage4SaveTime = datetime('now');
save(stage4MatFile);
writetable(resultTableStage4, stage4ParameterCsvFile);
writetable(corrTableStage4, stage4CorrelationCsvFile, 'WriteRowNames', true);

fprintf('\nSaved Stage 4 workspace and technical data:\n');
fprintf('  %s\n', stage4MatFile);
fprintf('  %s\n', stage4ParameterCsvFile);
fprintf('  %s\n', stage4CorrelationCsvFile);

if useDiary
    diary off;
end
