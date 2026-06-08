function referenceID = kalman_identify_reference_model(projectPaths, settings, referenceExperiments)
%KALMAN_IDENTIFY_REFERENCE_MODEL Identify or load the even-run reference model.

referenceSignature = kalman_reference_signature(settings);

if ~settings.forceReferenceIdentification && isfile(settings.referenceResultFile)
    cached = load(settings.referenceResultFile);
    if isfield(cached, 'referenceID') && isfield(cached.referenceID, 'signature') ...
            && isequaln(cached.referenceID.signature, referenceSignature)
        fprintf('Loaded cached Kalman reference model: %s\n', settings.referenceResultFile);
        referenceID = cached.referenceID;
        return;
    end
end

modelPath = which('greybox_id_full_stribeck_model');
if isempty(modelPath)
    error('Could not find greybox_id_full_stribeck_model.m on the MATLAB path.');
end

fprintf('\nIdentifying Stribeck reference model from even validation runs...\n');
fprintf('Using model file: %s\n', modelPath);

[zReference, initialStateValues, experimentNames] = kalman_make_iddata(referenceExperiments, settings);
passive = kalman_load_passive_stribeck_parameters(projectPaths.passiveLinkStribeckResult);

parameterNames = {'p_a', 'p_b1', 'p_c1', 'p_g1', 'p_u', 'p_0', ...
    'p_b2', 'p_g2', 'p_c2', 'p_sdelta2', 'v_s2', 'eps_v1', 'eps_v2'};

initialParameterValues = [
    14;
    670;
    0;
    97;
    4000;
    0;
    passive.p_b2;
    passive.p_g2;
    passive.p_c2;
    passive.p_sdelta2;
    passive.v_s2;
    0.05;
    passive.eps_v2
];

[initialParameterValues, initialParameterSource] = kalman_reference_initial_parameters( ...
    projectPaths, parameterNames, initialParameterValues);
parameters = num2cell(initialParameterValues);

minimumValues = [5, 100, 0, 5, 500, -200, 0, 0, 0, 0, 0.02, 0.001, 0.001];
maximumValues = [80, 4000, 800, 500, 30000, 200, Inf, Inf, Inf, Inf, 20, 1, 1];
fixedStage1 = [false, false, true, false, false, false, true, true, true, true, true, true, true];
fixedStage2 = [false, false, false, false, false, false, true, true, true, true, true, true, true];

initialStates = {
    initialStateValues(1, :);
    initialStateValues(2, :);
    initialStateValues(3, :);
    initialStateValues(4, :)
};

model0 = idnlgrey('greybox_id_full_stribeck_model', [2 1 4], parameters, initialStates, 0);
for k = 1:numel(parameterNames)
    model0.Parameters(k).Name = parameterNames{k};
    model0.Parameters(k).Minimum = minimumValues(k);
    model0.Parameters(k).Maximum = maximumValues(k);
    model0.Parameters(k).Fixed = fixedStage1(k);
end

model0.InputName = {'u'};
model0.InputUnit = {'V'};
model0.OutputName = {'theta_1', 'theta_2'};
model0.OutputUnit = {'rad', 'rad'};
model0.TimeUnit = 's';

model0.InitialStates(1).Name = 'theta_1';
model0.InitialStates(1).Unit = 'rad';
model0.InitialStates(2).Name = 'theta_2';
model0.InitialStates(2).Unit = 'rad';
model0.InitialStates(3).Name = 'theta_1_dot';
model0.InitialStates(3).Unit = 'rad/s';
model0.InitialStates(4).Name = 'theta_2_dot';
model0.InitialStates(4).Unit = 'rad/s';

nExp = numel(referenceExperiments);
model0 = setinit(model0, 'Fixed', {
    true(1, nExp);
    true(1, nExp);
    true(1, nExp);
    true(1, nExp)
});

optStage1 = nlgreyestOptions;
optStage1.Display = settings.referenceDisplay;
optStage1.EstimateCovariance = settings.estimateReferenceCovariance;
optStage1.SearchMethod = 'lm';
optStage1.SearchOptions.MaxIterations = settings.referenceStage1MaxIterations;
optStage1.OutputWeight = settings.referenceOutputWeight;

optStage2 = nlgreyestOptions;
optStage2.Display = settings.referenceDisplay;
optStage2.EstimateCovariance = settings.estimateReferenceCovariance;
optStage2.SearchMethod = 'lm';
optStage2.SearchOptions.MaxIterations = settings.referenceStage2MaxIterations;
optStage2.OutputWeight = settings.referenceOutputWeight;

fprintf('\nReference ID stage 1: active parameters with p_c1 fixed.\n');
modelStage1 = nlgreyest(zReference, model0, optStage1);
modelStage1.Name = 'Kalman reference model stage 1';

fprintf('\nReference ID stage 2: active parameters with p_c1 unlocked.\n');
modelStage2Start = modelStage1;
modelStage2Start.Parameters(3).Value = 10;
for k = 1:numel(parameterNames)
    modelStage2Start.Parameters(k).Fixed = fixedStage2(k);
end

referenceModel = nlgreyest(zReference, modelStage2Start, optStage2);
referenceModel.Name = 'Kalman tuning Stribeck reference model';

metadata = struct();
metadata.modelFile = modelPath;
metadata.passiveResultFile = projectPaths.passiveLinkStribeckResult;
metadata.experimentNames = experimentNames;
metadata.parameterNames = parameterNames;
metadata.parameters = kalman_model_parameter_vector(referenceModel);
metadata.initialParameterSource = initialParameterSource;
metadata.stage1Loss = modelStage1.Report.Fit.LossFcn;
metadata.stage2Loss = referenceModel.Report.Fit.LossFcn;

referenceID = struct();
referenceID.model = referenceModel;
referenceID.modelStage1 = modelStage1;
referenceID.zReference = zReference;
referenceID.signature = referenceSignature;
referenceID.metadata = metadata;

save(settings.referenceResultFile, 'referenceID');
fprintf('Saved Kalman reference model: %s\n', settings.referenceResultFile);
end
