function lin = linearize_rotpendulum(settings)
%LINEARIZE_ROTPENDULUM Linearize the selected no-p0 down-line pendulum model.
%
% The state is the measured/saved coordinate:
%   x = [theta1_meas; theta2_meas; theta1_dot_meas; theta2_dot_meas]
%
% If no x0 is provided, this function uses the measured-coordinate state
% corresponding to the physical all-down equilibrium:
%   theta1_phys = pi, theta2_phys = 0, velocities = 0.

if nargin < 1 || isempty(settings)
    settings = struct();
end
if isfield(settings, 'linearization')
    settings = settings.linearization;
end

scriptFolder = fileparts(mfilename('fullpath'));
addpath(scriptFolder, '-begin');

if isfield(settings, 'p')
    p = settings.p(:);
else
    p = load_parameters();
end

x0 = settingOrDefault(settings, 'x0', []);
if isempty(x0)
    x0 = measuredAllDownEquilibrium(p);
end

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

% With the selected no-p0 model, the physical equilibrium input at the
% all-down/down-line equilibrium is zero.
u0 = settingOrDefault(settings, 'u0', 0);

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
sys_lin.StateName = {'theta1_meas', 'theta2_meas', 'theta1_dot_meas', 'theta2_dot_meas'};
sys_lin.InputName = {'u'};
sys_lin.OutputName = {'theta1_meas', 'theta2_meas'};

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
lin.theta2_offset = pi - p(13)*p(15) - p(14);

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
    saveData.f0 = lin.f0;
    saveData.y0 = lin.y0;
    saveData.theta2_offset = lin.theta2_offset;
    saveData.settings = lin.settings;
    save(outputFile, '-struct', 'saveData');
end
end

function x_eq = measuredAllDownEquilibrium(p)
theta_scale        = max(p(13), 1e-6);
theta1_offset      = p(14);
theta_abs_down_raw = p(15);
theta2_offset      = pi - theta_scale*theta_abs_down_raw - theta1_offset;

theta1_meas_eq = (pi - theta1_offset) / theta_scale;
theta2_meas_eq = (0  - theta2_offset) / theta_scale;

x_eq = [theta1_meas_eq; theta2_meas_eq; 0; 0];
end

function value = settingOrDefault(settings, fieldName, defaultValue)
if isfield(settings, fieldName) && ~isempty(settings.(fieldName))
    value = settings.(fieldName);
else
    value = defaultValue;
end
end
