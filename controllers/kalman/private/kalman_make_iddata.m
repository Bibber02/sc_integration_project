function [z, initialStates, experimentNames] = kalman_make_iddata(experiments, settings)
%KALMAN_MAKE_IDDATA Convert experiment structs to multi-experiment iddata.

z = [];
nExp = numel(experiments);
initialStates = zeros(4, nExp);
experimentNames = cell(nExp, 1);

for k = 1:nExp
    y = experiments(k).y;
    u = experiments(k).u;

    zExp = iddata(y, u, settings.TsData);
    zExp.Name = experiments(k).name;
    zExp.InputName = {'u'};
    zExp.InputUnit = {'V'};
    zExp.OutputName = {'theta_1', 'theta_2'};
    zExp.OutputUnit = {'rad', 'rad'};
    zExp.TimeUnit = 's';

    if isempty(z)
        z = zExp;
    else
        z = merge(z, zExp);
    end

    initialStates(:, k) = [y(1, 1); y(1, 2); 0; 0];
    experimentNames{k} = experiments(k).name;
end
end

