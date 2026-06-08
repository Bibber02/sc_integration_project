function lin = linearize_rotpendulum(settings)
%LINEARIZE_ROTPENDULUM Linearize the reduced nonlinear pendulum model.

if nargin < 1 || isempty(settings)
    settings = defaultSettings();
else
    settings = normalizeSettings(settings);
end

scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = scriptFolder;
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(projectRoot);
scip.setupPath;
projectPaths = scip.paths;

run(fullfile(projectPaths.model, 'load_parameters.m'));

linSettings = settings.linearization;
x0 = linSettings.x0(:);
Ts = linSettings.Ts;
hx = linSettings.hx;
hu = linSettings.hu;
u0 = computeEquilibriumInput(p, linSettings);

[f0, y0] = nonlinearPlant(x0, u0, p);

disp('Operating point check:');
disp('x0 = ');
disp(x0);

disp('u0 = ');
disp(u0);

disp('f(x0,u0) = ');
disp(f0);

disp('y0 = ');
disp(y0);

if norm(f0) > 1e-6
    warning('The selected operating point is not an exact equilibrium. The linear model will have an affine offset.');
end

nx = numel(x0);
nu = 1;
ny = numel(y0);

A = zeros(nx,nx);
B = zeros(nx,nu);
C = zeros(ny,nx);
D = zeros(ny,nu);

for i = 1:nx
    dx = zeros(nx,1);
    dx(i) = hx;

    [f_plus,  y_plus]  = nonlinearPlant(x0 + dx, u0, p);
    [f_minus, y_minus] = nonlinearPlant(x0 - dx, u0, p);

    A(:,i) = (f_plus - f_minus)/(2*hx);
    C(:,i) = (y_plus - y_minus)/(2*hx);
end

[f_plus,  y_plus]  = nonlinearPlant(x0, u0 + hu, p);
[f_minus, y_minus] = nonlinearPlant(x0, u0 - hu, p);

B(:,1) = (f_plus - f_minus)/(2*hu);
D(:,1) = (y_plus - y_minus)/(2*hu);

sys_lin = ss(A,B,C,D);
sys_lin.StateName = {'theta1','theta2','theta1_dot','theta2_dot'};
sys_lin.InputName = {'u'};
sys_lin.OutputName = {'theta1','theta2'};

sys_disc = c2d(sys_lin, Ts, 'zoh');

Ad = sys_disc.A;
Bd = sys_disc.B;
Cd = sys_disc.C;
Dd = sys_disc.D;

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
lin.outputFile = linSettings.outputFile;

disp('Continuous-time model:');
disp('A = ');
disp(A);

disp('B = ');
disp(B);

disp('C = ');
disp(C);

disp('D = ');
disp(D);

disp('Eigenvalues of continuous-time A:');
disp(eig(A));

disp(' ');
disp('Discrete-time model:');
disp('Ad = ');
disp(Ad);

disp('Bd = ');
disp(Bd);

disp('Cd = ');
disp(Cd);

disp('Dd = ');
disp(Dd);

disp('Eigenvalues of discrete-time Ad:');
disp(eig(Ad));

if linSettings.saveOutput
    outputFolder = fileparts(linSettings.outputFile);
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    outputFile = linSettings.outputFile;
    save(outputFile, ...
         'A','B','C','D', ...
         'Ad','Bd','Cd','Dd', ...
         'sys_lin','sys_disc', ...
         'Ts','x0','u0','p','settings');

    fprintf('Saved continuous and discrete linearized models to %s\n', outputFile);
end
end

function settings = defaultSettings()
scriptFolder = fileparts(mfilename('fullpath'));
projectRoot = scriptFolder;
while ~isfolder(fullfile(projectRoot, '+scip')) && ~strcmp(projectRoot, fileparts(projectRoot))
    projectRoot = fileparts(projectRoot);
end
addpath(projectRoot);
scip.setupPath;
projectPaths = scip.paths;

settings.linearization.x0 = [pi; pi; 0; 0];
settings.linearization.u0Mode = 'cancel_input_bias';
settings.linearization.Ts = 0.001;
settings.linearization.hx = 1e-6;
settings.linearization.hu = 1e-6;
settings.linearization.saveOutput = true;
settings.linearization.outputFile = fullfile(projectPaths.model, 'linearized_plant.mat');
end

function settings = normalizeSettings(settings)
if isfield(settings, 'linearization')
    settings = mergeWithDefaults(settings);
    return;
end

wrapper.linearization = struct();

if isfield(settings, 'x0')
    wrapper.linearization.x0 = settings.x0;
end

if isfield(settings, 'sampleTime')
    wrapper.linearization.Ts = settings.sampleTime;
elseif isfield(settings, 'Ts')
    wrapper.linearization.Ts = settings.Ts;
end

if isfield(settings, 'u0')
    wrapper.linearization.u0 = settings.u0;
end

if isfield(settings, 'u0Mode')
    wrapper.linearization.u0Mode = settings.u0Mode;
end

if isfield(settings, 'hx')
    wrapper.linearization.hx = settings.hx;
end

if isfield(settings, 'hu')
    wrapper.linearization.hu = settings.hu;
end

if isfield(settings, 'saveOutput')
    wrapper.linearization.saveOutput = settings.saveOutput;
end

if isfield(settings, 'outputFile')
    wrapper.linearization.outputFile = settings.outputFile;
elseif isfield(settings, 'linearizedPlantFile')
    wrapper.linearization.outputFile = settings.linearizedPlantFile;
end

settings = mergeWithDefaults(wrapper);
end

function settings = mergeWithDefaults(settings)
defaults = defaultSettings();

if ~isfield(settings, 'linearization')
    settings.linearization = struct();
end

fields = fieldnames(defaults.linearization);
for k = 1:numel(fields)
    fieldName = fields{k};
    if ~isfield(settings.linearization, fieldName) || isempty(settings.linearization.(fieldName))
        settings.linearization.(fieldName) = defaults.linearization.(fieldName);
    end
end
end

function u0 = computeEquilibriumInput(p, linSettings)
if isfield(linSettings, 'u0')
    u0 = linSettings.u0;
    return;
end

switch linSettings.u0Mode
    case 'cancel_input_bias'
        u0 = -p(6)/p(5);
    otherwise
        error('Unknown input equilibrium mode: %s', linSettings.u0Mode);
end
end

function [xdot, y] = nonlinearPlant(x, u, p)
%#codegen
theta1  = x(1);
theta2  = x(2);
dtheta1 = x(3);
dtheta2 = x(4);

pa       = p(1);
pb1      = p(2);
pc1      = p(3);
pg1      = p(4);
pu       = p(5);
p0       = p(6);
pb2      = p(7);
pg2      = p(8);
pc2      = p(9);
psdelta2 = p(10);
vs2      = p(11);
epsv1    = p(12);
epsv2    = p(13);
pc       = p(14);

c2 = cos(theta2);
s2 = sin(theta2);

M11 = pa + 1 + 2*pc*c2;
M12 = 1 + pc*c2;
M21 = M12;
M22 = 1;

Cqdot1 = -pc*s2*(2*dtheta1*dtheta2 + dtheta2^2);
Cqdot2 =  pc*s2*dtheta1^2;

F1 = pb1*dtheta1 + pc1*tanh(dtheta1/epsv1);

stribeck2 = pc2 + psdelta2*exp(-(dtheta2/vs2)^2);
F2 = pb2*dtheta2 + stribeck2*tanh(dtheta2/epsv2);

G1 = -pg1*sin(theta1) - pg2*sin(theta1 + theta2);
G2 = -pg2*sin(theta1 + theta2);

I1 = pu * u + p0;
I2 = 0;

rhs1 = I1 - Cqdot1 - F1 - G1;
rhs2 = I2 - Cqdot2 - F2 - G2;

detM = M11*M22 - M12*M21;

ddtheta1 = ( M22*rhs1 - M12*rhs2) / detM;
ddtheta2 = (-M21*rhs1 + M11*rhs2) / detM;

xdot = zeros(4,1);
xdot(1) = dtheta1;
xdot(2) = dtheta2;
xdot(3) = ddtheta1;
xdot(4) = ddtheta2;

y = zeros(2,1);
y(1) = theta1;
y(2) = theta2;
end
