function lin = linearize_rotpendulum(settings)
%LINEARIZE_ROTPENDULUM Linearize the reduced nonlinear pendulum model.

if nargin < 1 || isempty(settings)
    settings = struct();
end
if isfield(settings, 'linearization')
    settings = settings.linearization;
end

scriptFolder = fileparts(mfilename('fullpath'));
addpath(scriptFolder, '-begin');

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

if isfield(settings, 'linearizedPlantFile') && isempty(settingOrDefault(settings, 'outputFile', []))
    outputFile = settings.linearizedPlantFile;
end

if isfield(settings, 'p')
    p = settings.p(:);
else
    p = load_parameters();
end

if isfield(settings, 'u0')
    u0 = settings.u0;
else
    u0 = -p(6) / p(5);
end

x0 = x0(:);
[f0, y0] = nonlinearPlant(x0, u0, p);

nx = numel(x0);
nu = 1;
ny = numel(y0);

A = zeros(nx, nx);
B = zeros(nx, nu);
C = zeros(ny, nx);
D = zeros(ny, nu);

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

sys_lin = ss(A, B, C, D);
sys_lin.StateName = {'theta1', 'theta2', 'theta1_dot', 'theta2_dot'};
sys_lin.InputName = {'u'};
sys_lin.OutputName = {'theta1', 'theta2'};

sys_disc = c2d(sys_lin, Ts, 'zoh');

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

if saveOutput
    outputFolder = fileparts(outputFile);
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    saveData = struct();
    saveData.A = lin.A;
    saveData.B = lin.B;
    saveData.C = lin.C;
    saveData.D = lin.D;
    saveData.Ad = lin.Ad;
    saveData.Bd = lin.Bd;
    saveData.Cd = lin.Cd;
    saveData.Dd = lin.Dd;
    saveData.sys_lin = lin.sys_lin;
    saveData.sys_disc = lin.sys_disc;
    saveData.Ts = lin.Ts;
    saveData.x0 = lin.x0;
    saveData.u0 = lin.u0;
    saveData.p = lin.p;
    saveData.settings = lin.settings;
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