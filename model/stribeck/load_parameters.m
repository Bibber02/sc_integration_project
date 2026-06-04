clear;
clc;

% ------------------------------------------------------------
% User settings
% ------------------------------------------------------------

matFile = 'identified_parameters.mat';   % <-- change to your actual .mat file name

% Measured distance from motor pivot to passive joint
l1 = 0.10;   % <-- replace with your real l1 value in meters

g = 9.81;

% ------------------------------------------------------------
% Load .mat file
% ------------------------------------------------------------

S = load(matFile);

% ------------------------------------------------------------
% Automatically find the 13x6 table inside the .mat file
% ------------------------------------------------------------

fieldNames = fieldnames(S);

T = [];
tableName = '';

for k = 1:numel(fieldNames)
    candidate = S.(fieldNames{k});

    if istable(candidate)
        if height(candidate) == 13 && width(candidate) >= 2
            T = candidate;
            tableName = fieldNames{k};
            break;
        end
    end
end

if isempty(T)
    error('No 13-row table with at least 2 columns was found in %s.', matFile);
end

fprintf('Loaded parameter table: %s\n', tableName);

% ------------------------------------------------------------
% Extract second column: parameter values
% ------------------------------------------------------------

values = T{:,2};

% Convert to numeric if needed
if iscell(values)
    values = cellfun(@double, values);
elseif isstring(values) || ischar(values)
    values = str2double(values);
end

values = double(values(:));

if numel(values) ~= 13
    error('Expected 13 parameter values, but found %d.', numel(values));
end

if any(isnan(values))
    error('At least one parameter value became NaN. Check the second column of the table.');
end

% ------------------------------------------------------------
% Expected row order
% ------------------------------------------------------------
% values(1)  = pa
% values(2)  = pb1
% values(3)  = pc1
% values(4)  = pg1
% values(5)  = pu
% values(6)  = p0
% values(7)  = pb2
% values(8)  = pg2
% values(9)  = pc2
% values(10) = psdelta2
% values(11) = vs2
% values(12) = epsv1
% values(13) = epsv2

pa       = values(1);
pb1      = values(2);
pc1      = values(3);
pg1      = values(4);
pu       = values(5);
%p0       = values(6);
p0       = 0;

pb2      = values(7);
pg2      = values(8);
pc2      = values(9);
psdelta2 = values(10);
vs2      = values(11);

epsv1    = values(12);
epsv2    = values(13);

% ------------------------------------------------------------
% Compute coupling parameter pc
% ------------------------------------------------------------
% In the reduced model:
%   pc = gamma / beta = (l1/g)*pg2
%
% Your ID table has 13 parameters, but the Simulink plant function
% also needs pc. Therefore we append pc as parameter 14.
% ------------------------------------------------------------

pc = (l1/g)*pg2;

% ------------------------------------------------------------
% Final parameter vector for Simulink
% ------------------------------------------------------------

p = [
    pa;
    pb1;
    pc1;
    pg1;
    pu;
    p0;
    pb2;
    pg2;
    pc2;
    psdelta2;
    vs2;
    epsv1;
    epsv2;
    pc
];

fprintf('Created Simulink parameter vector p with size %dx%d.\n', size(p,1), size(p,2));

% Optional: show result
disp(array2table(p, ...
    'VariableNames', {'Value'}, ...
    'RowNames', {'pa','pb1','pc1','pg1','pu','p0','pb2','pg2','pc2','psdelta2','vs2','epsv1','epsv2','pc'}));