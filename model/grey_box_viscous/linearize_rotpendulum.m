function lin = linearize_rotpendulum(settings)
%LINEARIZE_ROTPENDULUM Linearize the reduced nonlinear viscous pendulum model.

%% Default settings
if nargin < 1 || isempty(settings)
    settings = struct();
end

if isfield(settings, 'linearization')
    settings = settings.linearization;
end

%% Paths
scriptFolder = fileparts(mfilename('fullpath'));
addpath(scriptFolder, '-begin');

%% Read settings
x0 = settingOrDefault(settings, 'x0', [pi; pi; 0; 0]);

Ts = settingOrDefault(settings, 'sampleTime', []);
if isempty(Ts)
    Ts = settingOrDefault(settings, 'Ts', 0.01);
end

hx = settingOrDefault(settings, 'hx', 1e-6);
hu = settingOrDefault(settings, 'hu', 1e-6);

saveOutput = settingOrDefault(settings, 'saveOutput', false);
outputFile = settingOrDefault(settings, 'outputFile', ...
    fullfile(scriptFolder, 'linearized_plant.mat'));

%% Load parameters
if isfield(settings, 'p')
    p = settings.p(:);
else
    p = load_parameters_viscous();
end

%% Operating point input
if isfield(settings, 'u0')
    u0 = settings.u0;
else
    u0 = -p(5) / p(4);
end

x0 = x0(:);

%% Check operating point
[f0, y0] = nonlinearPlant(x0, u0, p);

nx = numel(x0);
nu = 1;
ny = numel(y0);

A = zeros(nx, nx);
B = zeros(nx, nu);
C = zeros(ny, nx);
D = zeros(ny, nu);

%% Linearize with central finite differences
for i = 1:nx
    dx = zeros(nx, 1);
    dx(i) = hx;

    [fPlus, yPlus] = nonlinearPlant(x0 + dx, u0, p);
    [fMinus, yMinus] = nonlinearPlant(x0 - dx, u0, p);

    A(:, i) = (fPlus - fMinus) / (2 * hx);
    C(:, i) = (yPlus - yMinus) / (2 * hx);
end

[fPlus, yPlus] = nonlinearPlant(x0, u0 + hu, p);
[fMinus, yMinus] = nonlinearPlant(x0, u0 - hu, p);

B(:, 1) = (fPlus - fMinus) / (2 * hu);
D(:, 1) = (yPlus - yMinus) / (2 * hu);

%% State-space models
sys_lin = ss(A, B, C, D);

sys_lin.StateName = {'theta1', 'theta2', 'theta1_dot', 'theta2_dot'};
sys_lin.InputName = {'u'};
sys_lin.OutputName = {'theta1', 'theta2'};

sys_disc = c2d(sys_lin, Ts, 'zoh');

%% Store output
lin = struct();

lin.A = A;
lin.B = B;
lin.C = C;
lin.D = D;

lin.Ad = sys_disc.A;
lin.Bd = sys_disc.B;
lin.Cd = sys_disc.C;
lin.Dd = sys_disc.D;

lin.sys_lin = sys_lin;
lin.sys_disc = sys_disc;

lin.Ts = sys_disc.Ts;
lin.x0 = x0;
lin.u0 = u0;
lin.p = p;

lin.f0 = f0;
lin.y0 = y0;

lin.settings = settings;
lin.outputFile = outputFile;

%% Optional save
if saveOutput
    outputFolder = fileparts(outputFile);

    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    saveData = lin;
    save(outputFile, '-struct', 'saveData');
end

end

function value = settingOrDefault(settings, fieldName, defaultValue)
if isfield(settings, fieldName) && ~isempty(settings.(fieldName))
    value = settings.(fieldName);
else
    value = defaultValue;
end
end
