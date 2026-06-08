if ~exist('sys_disc', 'var')
    error('calc_lqr requires sys_disc in the workspace. Run controllers/lqr/setup_lqr.m first.');
end

if exist('Q_lqr', 'var') ~= 1 || exist('R_lqr', 'var') ~= 1
    if exist('lqrSettings', 'var') == 1 && isstruct(lqrSettings)
        if isfield(lqrSettings, 'lqr')
            Q_lqr = lqrSettings.lqr.Q_lqr;
            R_lqr = lqrSettings.lqr.R_lqr;
        else
            Q_lqr = lqrSettings.Q_lqr;
            R_lqr = lqrSettings.R_lqr;
        end
    elseif exist('settings', 'var') == 1 && isstruct(settings)
        Q_lqr = settings.Q_lqr;
        R_lqr = settings.R_lqr;
    else
        error('calc_lqr requires Q_lqr and R_lqr in the workspace. Run controllers/lqr/setup_lqr.m first.');
    end
end

% Extract discrete matrices from the LTI system
Ad = sys_disc.A;
Bd = sys_disc.B;
Cd = sys_disc.C;
Dd = sys_disc.D;
Ts = sys_disc.Ts;

% Discrete-time LQR gain for model coordinates:
%   u_model_dev = -K_lqr * x_dev.
% The real hardware command may need an extra sign conversion; setup_lqr.m
% computes K_lqr_command for the Simulink Gain block.
K_lqr = dlqr(Ad, Bd, Q_lqr, R_lqr);
closedLoopPoles = eig(Ad - Bd*K_lqr);

disp('K_lqr = ');
disp(K_lqr);

disp('Closed-loop poles:');
disp(closedLoopPoles);

fs = 1/Ts;
filterOrder = 1;
fc = 10;

if exist('hiddenDefaults', 'var') == 1 && isstruct(hiddenDefaults)
    if isfield(hiddenDefaults, 'filterOrder')
        filterOrder = hiddenDefaults.filterOrder;
    end
    if isfield(hiddenDefaults, 'filterCutoffHz')
        fc = hiddenDefaults.filterCutoffHz;
    end
elseif exist('lqrSettings', 'var') == 1 && isstruct(lqrSettings) && isfield(lqrSettings, 'filter')
    if isfield(lqrSettings.filter, 'order')
        filterOrder = lqrSettings.filter.order;
    end
    if isfield(lqrSettings.filter, 'cutoffHz')
        fc = lqrSettings.filter.cutoffHz;
    end
end

wn = fc / (fs/2);
if wn >= 1
    error('Filter cutoff %.3f Hz must be below Nyquist frequency %.3f Hz.', fc, fs/2);
end
[b, a] = butter(filterOrder, wn);
