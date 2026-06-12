function p = load_parameters()
%LOAD_PARAMETERS Load the identified full-system parameters.

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptFolder));

matFile = fullfile(projectRoot, ...
    'system_identification', 'full_system', 'grey_box', 'stribeck', ...
    'results_with_p0', 'full_system_id_new_data_locked_passive_result.mat');


l1 = 0.10;
g = 9.81;

S = load(matFile, 'resultTable');
values = double(S.resultTable{:, 2});

pa       = values(1);
pb1      = values(2);
pc1      = values(3);
pg1      = values(4);
pu       = values(5);
p0       = 0;
pb2      = values(7);
pg2      = values(8);
pc2      = values(9);
psdelta2 = values(10);
vs2      = values(11);
epsv1    = values(12);
epsv2    = values(13);
pc       = (l1 / g) * pg2;

p = [
    pa
    pb1
    pc1
    pg1
    pu
    p0
    pb2
    pg2
    pc2
    psdelta2
    vs2
    epsv1
    epsv2
    pc
];

fprintf('Loaded parameter vector p from %s\n', matFile);
end