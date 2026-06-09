%% Single-file LQR setup
% Edit this block, then run this file from anywhere.

sampleTime = 0.001;

% Angle convention:
%   all-down: [pi; 0]
%   down-up:  [pi; pi]
%   all-up:   [0; 0]
x0 = [pi; pi; 0; 0];

Q_lqr = diag([5 1 0.1 0.1]);
R_lqr = 5;

filterOrder = 1;
filterCutoffHz = 10;

% Parameter loading / linearization settings.
l1 = 0.10;
g = 9.81;
hx = 1e-6;
hu = 1e-6;

% Hardware calibration constants used by the Simulink I/O subsystem.
daoutoffs = 0.00;
daoutgain = -6;
adinoffs = -[1.070535559513546 1.190099998448171];
adingain = [1.200374077783257 1.222725978474668];
adinoffs = [adinoffs 0 0 0 0 0];
adingain = [adingain 1 1 1 1 1];
sensorChannels = [6 7];
hardwareInitEnabled = true;

% Kalman filter settings. These can be replaced by tuned values later.
Q_kalman = diag([0.1 0.1 1 0.1]);
R_kalman = diag([7.3339511547e-6 1.13475782135e-5]);
P0_kalman = 10 * eye(4);
x0_kalman = zeros(4, 1);

openModelAfterSetup = true;

%% Setup logic

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptFolder));
templateSupportFolder = fullfile(projectRoot, 'template', 'rotating-pendulum');
parameterResultFile = fullfile(projectRoot, ...
    'system_identification', 'full_system', 'grey_box', 'stribeck', ...
    'results_with_p0', 'full_system_id_new_data_locked_passive_result.mat');
linearizedPlantFile = fullfile(scriptFolder, 'lqr_linearized_plant.mat');
setupResultFile = fullfile(scriptFolder, 'lqr_setup_result.mat');
simulinkModel = fullfile(scriptFolder, 'rotating-pendulum', 'rotpentemplate.slx');

cd(projectRoot);
addpath(scriptFolder, '-begin');

settings = struct();
settings.sampleTime = sampleTime;
settings.x0 = x0;
settings.Q_lqr = Q_lqr;
settings.R_lqr = R_lqr;
settings.filterOrder = filterOrder;
settings.filterCutoffHz = filterCutoffHz;
settings.l1 = l1;
settings.g = g;
settings.hx = hx;
settings.hu = hu;
settings.parameterResultFile = parameterResultFile;
settings.linearizedPlantFile = linearizedPlantFile;
settings.setupResultFile = setupResultFile;
settings.simulinkModel = simulinkModel;
settings.hardwareInitEnabled = hardwareInitEnabled;

lqrSettings = settings;
hiddenDefaults = struct();
hiddenDefaults.filterOrder = filterOrder;
hiddenDefaults.filterCutoffHz = filterCutoffHz;

p = loadPlantParameters(parameterResultFile, l1, g);
u0 = -p(6) / p(5);

[A, B, C, D, y0, f0] = linearizePlant(x0, u0, p, hx, hu);

sys_lin = ss(A, B, C, D);
sys_lin.StateName = {'theta1', 'theta2', 'theta1_dot', 'theta2_dot'};
sys_lin.InputName = {'u'};
sys_lin.OutputName = {'theta1', 'theta2'};

sys_disc = c2d(sys_lin, sampleTime, 'zoh');
Ad = sys_disc.A;
Bd = sys_disc.B;
Cd = sys_disc.C;
Dd = sys_disc.D;
Ts = sys_disc.Ts;
h = Ts;
runHardwareInit = false; %#ok<NASGU>

run(fullfile(scriptFolder, 'calc_lqr.m'));

lin = struct();
lin.A = A;
lin.B = B;
lin.C = C;
lin.D = D;
lin.Ad = Ad;
lin.Bd = Bd;
lin.Cd = Cd;
lin.Dd = Dd;
lin.sys_lin = sys_lin;
lin.sys_disc = sys_disc;
lin.Ts = Ts;
lin.x0 = x0;
lin.u0 = u0;
lin.p = p;
lin.f0 = f0;
lin.y0 = y0;
lin.settings = settings;
lin.outputFile = linearizedPlantFile;

save(linearizedPlantFile, ...
    'A', 'B', 'C', 'D', 'Ad', 'Bd', 'Cd', 'Dd', ...
    'sys_lin', 'sys_disc', 'Ts', 'x0', 'u0', 'p', 'settings');

save(setupResultFile, ...
    'settings', 'hiddenDefaults', 'lin', ...
    'A', 'B', 'C', 'D', 'Ad', 'Bd', 'Cd', 'Dd', ...
    'sys_lin', 'sys_disc', 'Ts', 'x0', 'y0', 'u0', 'p', ...
    'Q_lqr', 'R_lqr', 'K_lqr', 'closedLoopPoles', 'b', 'a', ...
    'Q_kalman', 'R_kalman', 'P0_kalman', 'x0_kalman', ...
    'h', 'daoutoffs', 'daoutgain', 'adinoffs', 'adingain', ...
    'sensorChannels', 'hardwareInitEnabled');

workspaceVars = {
    'settings'
    'lqrSettings'
    'hiddenDefaults'
    'lin'
    'A'
    'B'
    'C'
    'D'
    'Ad'
    'Bd'
    'Cd'
    'Dd'
    'sys_lin'
    'sys_disc'
    'Ts'
    'h'
    'x0'
    'y0'
    'u0'
    'p'
    'Q_lqr'
    'R_lqr'
    'K_lqr'
    'closedLoopPoles'
    'b'
    'a'
    'Q_kalman'
    'R_kalman'
    'P0_kalman'
    'x0_kalman'
    'daoutoffs'
    'daoutgain'
    'adinoffs'
    'adingain'
    'sensorChannels'
    'hardwareInitEnabled'
    'runHardwareInit'
};

for k = 1:numel(workspaceVars)
    assignin('base', workspaceVars{k}, eval(workspaceVars{k}));
end

load_system(simulinkModel);
[~, modelName] = fileparts(simulinkModel);
addpath(templateSupportFolder, '-begin');
configureSimulinkModel(modelName, sensorChannels);

if openModelAfterSetup && usejava('desktop')
    open_system(modelName);
else
    close_system(modelName, 0);
end

runHardwareInit = hardwareInitEnabled;
assignin('base', 'runHardwareInit', runHardwareInit);

fprintf('\nLQR setup complete.\n');
fprintf('Project root:      %s\n', projectRoot);
fprintf('Sample time:       %.6g s\n', Ts);
fprintf('Initial state x0:  [%g; %g; %g; %g]\n', x0(1), x0(2), x0(3), x0(4));
fprintf('Linearized plant:  %s\n', linearizedPlantFile);
fprintf('Setup result:      %s\n', setupResultFile);
fprintf('Simulink model:    %s\n', simulinkModel);

function p = loadPlantParameters(parameterResultFile, l1, g)
S = load(parameterResultFile, 'resultTable');
values = S.resultTable{:, 2};

if iscell(values)
    values = cellfun(@double, values);
elseif isstring(values) || ischar(values)
    values = str2double(values);
end

values = double(values(:));

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
end

function [A, B, C, D, y0, f0] = linearizePlant(x0, u0, p, hx, hu)
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
end

function [xdot, y] = nonlinearPlant(x, u, p)
theta1 = x(1);
theta2 = x(2);
dtheta1 = x(3);
dtheta2 = x(4);

pa = p(1);
pb1 = p(2);
pc1 = p(3);
pg1 = p(4);
pu = p(5);
p0 = p(6);
pb2 = p(7);
pg2 = p(8);
pc2 = p(9);
psdelta2 = p(10);
vs2 = p(11);
epsv1 = p(12);
epsv2 = p(13);
pc = p(14);

c2 = cos(theta2);
s2 = sin(theta2);

M11 = pa + 1 + 2 * pc * c2;
M12 = 1 + pc * c2;
M21 = M12;
M22 = 1;

Cqdot1 = -pc * s2 * (2 * dtheta1 * dtheta2 + dtheta2^2);
Cqdot2 = pc * s2 * dtheta1^2;

F1 = pb1 * dtheta1 + pc1 * tanh(dtheta1 / epsv1);

stribeck2 = pc2 + psdelta2 * exp(-(dtheta2 / vs2)^2);
F2 = pb2 * dtheta2 + stribeck2 * tanh(dtheta2 / epsv2);

G1 = -pg1 * sin(theta1) - pg2 * sin(theta1 + theta2);
G2 = -pg2 * sin(theta1 + theta2);

I1 = pu * u + p0;
I2 = 0;

rhs1 = I1 - Cqdot1 - F1 - G1;
rhs2 = I2 - Cqdot2 - F2 - G2;

detM = M11 * M22 - M12 * M21;

ddtheta1 = (M22 * rhs1 - M12 * rhs2) / detM;
ddtheta2 = (-M21 * rhs1 + M11 * rhs2) / detM;

xdot = zeros(4, 1);
xdot(1) = dtheta1;
xdot(2) = dtheta2;
xdot(3) = ddtheta1;
xdot(4) = ddtheta2;

y = zeros(2, 1);
y(1) = theta1;
y(2) = theta2;
end

function configureSimulinkModel(modelName, sensorChannels)
set_param([modelName '/Gain'], 'Gain', '-K_lqr');
set_param([modelName '/Constant'], 'Value', 'y0');

set_param([modelName '/Kalman Filter'], ...
    'ModelSourceVariable', 'sys_disc', ...
    'Q', 'Q_kalman', ...
    'R', 'R_kalman', ...
    'P0', 'P0_kalman', ...
    'X0', 'x0_kalman');

set_param([modelName '/Subsystem/RT input'], ...
    'Ts', 'h', ...
    'channels', mat2str(sensorChannels), ...
    'MaskInitialization', hardwareMaskInitialization());
end

function code = hardwareMaskInitialization()
code = sprintf([ ...
    'Ts = evalin(''base'', ''Ts'');\n' ...
    'h = evalin(''base'', ''h'');\n' ...
    'daoutoffs = evalin(''base'', ''daoutoffs'');\n' ...
    'daoutgain = evalin(''base'', ''daoutgain'');\n' ...
    'adinoffs = evalin(''base'', ''adinoffs'');\n' ...
    'adingain = evalin(''base'', ''adingain'');\n' ...
    'runHardwareInit = evalin(''base'', ''exist(''''runHardwareInit'''', ''''var'''') && runHardwareInit'');\n' ...
    'if runHardwareInit\n' ...
    'fugihandle = fugiboard(''Open'', ''Pendulum1'');\n' ...
    'fugihandle.WatchdogTimeout = 0.5;\n' ...
    'fugiboard(''SetParams'', fugihandle);\n' ...
    'fugiboard(''Write'', fugihandle, 0, 0, 0, 0);\n' ...
    'fugiboard(''Write'', fugihandle, 5, 1, 0, 0);\n' ...
    'data = fugiboard(''Read'', fugihandle);\n' ...
    'fpgaModel = bitshift(data(1), -5);\n' ...
    'fpgaVersion = bitand(data(1), 31);\n' ...
    'disp(sprintf(''FPGA setup %%d, version %%d'', fpgaModel, fpgaVersion));\n' ...
    'fugiboard(''Write'', fugihandle, 0, 1, 0, 0);\n' ...
    'pause(0.1);\n' ...
    'else\n' ...
    'fugihandle = [];\n' ...
    'end\n']);
end
