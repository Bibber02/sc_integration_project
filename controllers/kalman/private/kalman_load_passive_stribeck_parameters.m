function passive = kalman_load_passive_stribeck_parameters(passiveResultMatFile)
%KALMAN_LOAD_PASSIVE_STRIBECK_PARAMETERS Load fixed passive Stribeck params.

if ~isfile(passiveResultMatFile)
    error('Passive-link Stribeck result file not found: %s', passiveResultMatFile);
end

S = load(passiveResultMatFile);
passiveModel = [];

if isfield(S, 'passiveID')
    candidates = {'mStrBest', 'mCoul', 'mVisc', 'model_viscous'};
    for k = 1:numel(candidates)
        if isfield(S.passiveID, candidates{k})
            passiveModel = S.passiveID.(candidates{k});
            break;
        end
    end
else
    candidates = {'mStrBest', 'mCoul', 'mVisc', 'model_viscous'};
    for k = 1:numel(candidates)
        if isfield(S, candidates{k})
            passiveModel = S.(candidates{k});
            break;
        end
    end
end

if isempty(passiveModel)
    error('Could not find a passive-link idnlgrey model in %s.', passiveResultMatFile);
end

passive = struct('p_b2', NaN, 'p_g2', NaN, 'p_c2', NaN, ...
    'p_sdelta2', NaN, 'v_s2', NaN, 'eps_v2', NaN);

for k = 1:numel(passiveModel.Parameters)
    name = char(passiveModel.Parameters(k).Name);
    value = passiveModel.Parameters(k).Value;
    if iscell(value)
        value = value{1};
    end
    value = double(value);

    switch name
        case 'p_b2'
            passive.p_b2 = value;
        case 'p_g2'
            passive.p_g2 = value;
        case 'p_c2'
            passive.p_c2 = value;
        case {'p_sdelta2', 'p_sdelta'}
            passive.p_sdelta2 = value;
        case {'v_s2', 'v_s'}
            passive.v_s2 = value;
        case {'eps_v2', 'eps_v'}
            passive.eps_v2 = value;
    end
end

values = [passive.p_b2, passive.p_g2, passive.p_c2, ...
    passive.p_sdelta2, passive.v_s2, passive.eps_v2];
if any(isnan(values))
    error('Passive-link Stribeck result is missing one or more required parameters.');
end
end

