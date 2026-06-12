function p = load_parameters_viscous()
%LOAD_PARAMETERS_VISCOUS Load the identified viscous full-system parameters.

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptFolder));

matFile = fullfile(projectRoot, ...
    'system_identification', 'full_system', 'grey_box', 'viscous', ...
    'full_system_id_viscous_result.mat');

l1 = 0.10;
g = 9.81;

S = load(matFile, 'resultTable');
resultTable = S.resultTable;

requiredNames = {'p_a', 'p_b1', 'p_g1', 'p_u', 'p_0', 'p_b2', 'p_g2'};
values = zeros(numel(requiredNames), 1);

parameterNames = string(resultTable.Parameter);
parameterKeys = lower(regexprep(parameterNames, '[^a-zA-Z0-9]', ''));

for k = 1:numel(requiredNames)
    requiredKey = lower(regexprep(requiredNames{k}, '[^a-zA-Z0-9]', ''));
    match = find(parameterKeys == requiredKey, 1);

    if isempty(match)
        error('Could not find parameter "%s" in %s.', requiredNames{k}, matFile);
    end

    values(k) = double(resultTable.Value(match));
end

pa  = values(1);
pb1 = values(2);
pg1 = values(3);
pu  = values(4);
p0  = 0;
pb2 = values(6);
pg2 = values(7);
pc  = (l1 / g) * pg2;

p = [
    pa
    pb1
    pg1
    pu
    p0
    pb2
    pg2
    pc
];

fprintf('Loaded viscous parameter vector p from %s\n', matFile);
end
