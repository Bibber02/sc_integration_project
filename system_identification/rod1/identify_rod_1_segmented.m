%% ============================================================
%  Rod-1 (and motor) parameter identification — multi-segment
%
%  Chops each long PRBS recording into short overlapping segments
%  and treats them as separate experiments in a merged iddata.
%  Initial states (positions AND velocities) are taken from the
%  data via Savitzky-Golay filtering and held FIXED — no extra
%  free parameters for nlgreyest beyond the physical ones.
%
%  Input sign flipped at the iddata-building step so that p_u
%  stays positive in the model.
%
%  Assumes the workspace contains expData with .train and .val
%  sub-structs (.t, .u, .theta_1, .theta_2) per experiment.
% ============================================================

%% ----- Sanity checks ----------------------------------------
if ~exist('expData', 'var')
    error('expData not found. Run the data-loading cells of identificatino_rod_1.mlx first.');
end
assert(numel(expData) >= 2, 'Expected at least two experiments.');

Ts = expData(1).train.t(2) - expData(1).train.t(1);
fprintf('Sampling time Ts = %g s\n', Ts);

%% ----- Segmentation settings --------------------------------
T_seg     = 10.0;     % [s] segment length
overlap   = 0.5;      % 50% overlap
T_stride  = T_seg * (1 - overlap);
minSegFrac = 0.8;     % keep segments at least 80% of T_seg

fprintf('Segment length: %.1f s, stride: %.1f s (overlap %.0f%%)\n', ...
        T_seg, T_stride, 100*overlap);

%% ----- Fixed (known) constants ------------------------------
l_1     = 0.100;
g_const = 9.81;

% Rod-2 from labmate's fit (low-velocity damping model)
p_b2    = 0.0215307;
p_g2    = 112.016;
p_c2    = 0.464082;
p_b_low = 0.244986;
v_b     = 0.514232;
eps_v2  = 0.02;

%% ----- Initial guesses for rod-1 parameters -----------------
% From hand-tuning with reduced-damping settings.
% Note: input sign is flipped at iddata level so p_u is POSITIVE here.

p_alpha_0 = 18.0;
p_g1_0    = 150.0;
p_b1_0    = 0.7;
p_c1_0    = 0.5;
eps_v1    = 0.02;
p_u_0     = 180.0;
p_0_0     = 0.0;

p0 = [
    p_alpha_0;          % 1
    p_g1_0;             % 2
    p_b1_0;             % 3
    p_c1_0;             % 4
    eps_v1;             % 5  (fixed)
    p_u_0;              % 6
    p_0_0;              % 7
    l_1;                % 8  (fixed)
    p_b2;               % 9  (fixed)
    p_g2;               % 10 (fixed)
    p_c2;               % 11 (fixed)
    p_b_low;            % 12 (fixed)
    v_b;                % 13 (fixed)
    eps_v2;             % 14 (fixed)
    g_const             % 15 (fixed)
];

paramNames = {'p_alpha';'p_g1';'p_b1';'p_c1';'eps_v1'; ...
              'p_u';'p_0';'l_1'; ...
              'p_b2';'p_g2';'p_c2';'p_b_low';'v_b';'eps_v2';'g'};

paramUnits = {'-';'1/s^2';'1/s';'rad/s^2';'rad/s'; ...
              '1/s^2';'1/s^2';'m'; ...
              '1/s';'1/s^2';'rad/s^2';'1/s';'rad/s';'rad/s';'m/s^2'};

parameterFixed = [
    false; false; false; false;   true; ...   % p_alpha, p_g1, p_b1, p_c1, eps_v1
    false; false;                  true; ...  % p_u, p_0, l_1
    true; true; true; true; true; true; true  % rod-2 + g
];

paramMin = [ 0.1; 0;   0;   0;   1e-4;   0;   -50; 0;     0; 0; 0; 0; 0.005; 1e-4; 0   ];
paramMax = [ 100; Inf; Inf; Inf; 0.05;   1000;  50; 1;   Inf; Inf; Inf; Inf; 1.0; 0.05; 20 ];

%% ----- Build segmented training iddata ----------------------
fprintf('\nBuilding training segments...\n');
[z_train, segInfoTrain] = buildSegmentedIddata( ...
    expData, 'train', T_seg, T_stride, minSegFrac, Ts);

fprintf('  Training: %d segments total\n', numel(segInfoTrain));

% Initial states for all training segments (fixed, from data)
x0_train = zeros(4, numel(segInfoTrain));
for k = 1:numel(segInfoTrain)
    x0_train(:,k) = segInfoTrain(k).x0;
end

%% ----- Build segmented validation iddata --------------------
fprintf('Building validation segments...\n');
[z_val, segInfoVal] = buildSegmentedIddata( ...
    expData, 'val', T_seg, T_stride, minSegFrac, Ts);
fprintf('  Validation: %d segments total\n', numel(segInfoVal));

x0_val = zeros(4, numel(segInfoVal));
for k = 1:numel(segInfoVal)
    x0_val(:,k) = segInfoVal(k).x0;
end

%% ----- Build the idnlgrey model -----------------------------
FileName = 'pendulum_rod_1_model';
Order    = [2 1 4];
Ts_model = 0;

Parameters = num2cell(p0(:));
model = idnlgrey(FileName, Order, Parameters, x0_train, Ts_model);

for j = 1:numel(paramNames)
    model.Parameters(j).Name    = paramNames{j};
    model.Parameters(j).Unit    = paramUnits{j};
    model.Parameters(j).Minimum = paramMin(j);
    model.Parameters(j).Maximum = paramMax(j);
    model.Parameters(j).Fixed   = parameterFixed(j);
end

% Hold ALL initial states fixed (Option A)
nExp = size(x0_train, 2);
model = setinit(model, 'Name', {'theta_1';'theta_2';'theta_1_dot';'theta_2_dot'});
model = setinit(model, 'Unit', {'rad';'rad';'rad/s';'rad/s'});
model = setinit(model, 'Minimum', repmat({-Inf(1,nExp)}, 4, 1));
model = setinit(model, 'Maximum', repmat({ Inf(1,nExp)}, 4, 1));
model = setinit(model, 'Fixed',   repmat({ true(1,nExp)}, 4, 1));

%% ----- Fit ---------------------------------------------------
opt = nlgreyestOptions;
opt.Display = 'on';
opt.SearchMethod = 'auto';
opt.SearchOptions.MaxIterations = 200;

fprintf('\nStarting nlgreyest on %d segments (%.0f s each, total %.0f s of data)...\n', ...
        nExp, T_seg, nExp*T_stride);
estimated_model = nlgreyest(z_train, model, opt);

%% ----- Report parameter estimates ---------------------------
fprintf('\nEstimated rod-1 parameters:\n');
fprintf('---------------------------------------------------------------\n');
fprintf('%-10s %-14s %-14s %-8s\n', 'Name', 'Value', 'StdDev', 'Status');
fprintf('---------------------------------------------------------------\n');

% Pre-compute covariance once
try
    covar = getcov(estimated_model);
catch
    covar = [];
end
freeMask = ~arrayfun(@(p) p.Fixed, estimated_model.Parameters);
freeIdx  = cumsum(freeMask);

for j = 1:numel(estimated_model.Parameters)
    P = estimated_model.Parameters(j);
    if P.Fixed
        sd_str = '   (fixed)';
        status = 'fixed';
    else
        if ~isempty(covar) && freeIdx(j) <= size(covar,1)
            sd_val = sqrt(covar(freeIdx(j), freeIdx(j)));
            sd_str = sprintf('%.4g', sd_val);
        else
            sd_str = 'n/a';
        end
        status = 'est';
    end
    fprintf('%-10s %-14.6g %-14s %-8s\n', P.Name, P.Value, sd_str, status);
end
fprintf('---------------------------------------------------------------\n');

%% ----- Sanity check: implied natural frequencies ------------
p_alpha_hat = estimated_model.Parameters(1).Value;
p_g1_hat    = estimated_model.Parameters(2).Value;
p_c_hat     = (l_1/g_const) * p_g2;
omega_1_sq  = p_g1_hat / (p_alpha_hat + 1 + 2*p_c_hat);
omega_2_sq  = p_g2;
fprintf('\nImplied small-signal natural frequencies (rad/s):\n');
fprintf('  omega_1 = %.3f  (period %.2f s)\n', sqrt(abs(omega_1_sq)), 2*pi/sqrt(abs(omega_1_sq)));
fprintf('  omega_2 = %.3f  (period %.2f s)\n', sqrt(omega_2_sq), 2*pi/sqrt(omega_2_sq));

%% ----- Validation -------------------------------------------
% Apply estimated parameters to a model with val-segment initial states
val_model = estimated_model;
nExpVal = size(x0_val, 2);
val_model.InitialStates = [];
val_model = idnlgrey(FileName, Order, ...
    arrayfun(@(p) p.Value, estimated_model.Parameters, 'UniformOutput', false), ...
    x0_val, Ts_model);
for j = 1:numel(paramNames)
    val_model.Parameters(j).Name  = estimated_model.Parameters(j).Name;
    val_model.Parameters(j).Unit  = estimated_model.Parameters(j).Unit;
    val_model.Parameters(j).Fixed = true;
end
val_model = setinit(val_model, 'Fixed', repmat({true(1,nExpVal)}, 4, 1));

%% ----- Comparison plots: per-experiment overlay -------------
% Recombine segments into full traces, simulate end-to-end with
% fitted parameters using ONE initial state per original experiment.
% This is the more meaningful visualisation.

p_hat = arrayfun(@(p) p.Value, estimated_model.Parameters);

plotFullTraceComparison(expData, 'train', p_hat, FileName);
plotFullTraceComparison(expData, 'val',   p_hat, FileName);

%% ============================================================
%  Helper functions
%% ============================================================

function [zMerged, segInfo] = buildSegmentedIddata(expData, fieldName, T_seg, T_stride, minSegFrac, Ts)
    % Chop each experiment's recording into overlapping segments,
    % flip input sign, and merge into a single multi-experiment iddata.

    segCells = {};
    segInfo  = [];

    for i = 1:numel(expData)
        d  = expData(i).(fieldName);
        t  = d.t(:);   t = t - t(1);
        u  = d.u(:);
        th1 = d.theta_1(:);
        th2 = d.theta_2(:);

        nSeg = floor((t(end) - T_seg) / T_stride) + 1;
        nSeg = max(nSeg, 1);

        for s = 1:nSeg
            t_start = (s-1) * T_stride;
            t_end   = t_start + T_seg;

            idxMask = (t >= t_start) & (t < t_end);
            if nnz(idxMask) < minSegFrac * (T_seg / Ts)
                continue;  % skip short tail
            end

            t_seg   = t(idxMask);
            u_seg   = u(idxMask);
            th1_seg = th1(idxMask);
            th2_seg = th2(idxMask);

            % Normalise time so each segment starts at 0
            t_seg = t_seg - t_seg(1);

            % Flip input sign so that p_u is POSITIVE in the model
            u_seg = -u_seg;

            % Build iddata for this segment
            z = iddata([th1_seg, th2_seg], u_seg, Ts);
            z.InputName  = {'u'};
            z.InputUnit  = {'-'};
            z.OutputName = {'theta_1','theta_2'};
            z.OutputUnit = {'rad','rad'};
            z.TimeUnit   = 's';
            z.ExperimentName = sprintf('%s_%s_seg%02d', ...
                                       expData(i).name, fieldName, s);

            % Initial state from data (Savitzky-Golay velocity estimate)
            x0 = estimateInitialStateSG(t_seg, th1_seg, th2_seg);

            segCells{end+1} = z; %#ok<AGROW>
            info.expIdx  = i;
            info.segIdx  = s;
            info.x0      = x0;
            info.t_start = t_start;
            info.name    = z.ExperimentName;
            segInfo = [segInfo, info]; %#ok<AGROW>
        end
    end

    if isempty(segCells)
        error('No segments produced. Check T_seg / data length.');
    end

    zMerged = merge(segCells{:});
end


function x0 = estimateInitialStateSG(t, th1, th2)
    % Initial position from sample 1; initial velocity from a
    % Savitzky-Golay (order 2, window 11) derivative estimate
    % evaluated at sample 1. Robust to potentiometer noise.

    N = numel(t);
    win = min(11, 2*floor((N-1)/2)+1);   % ensure odd, <= N
    win = max(win, 3);

    if win < 5
        % Fallback: short polyfit on first few samples
        Nfit = min(5, N);
        c1 = polyfit(t(1:Nfit), th1(1:Nfit), 1);
        c2 = polyfit(t(1:Nfit), th2(1:Nfit), 1);
        v1 = c1(1);
        v2 = c2(1);
    else
        [~, gSG] = sgolay(2, win);  % gSG(:,2) is the 1st-derivative kernel
        Ts = t(2) - t(1);
        dth1 = conv(th1, factorial(1) * gSG(:,2) / (-Ts)^1, 'same');
        dth2 = conv(th2, factorial(1) * gSG(:,2) / (-Ts)^1, 'same');
        % Use the value at the middle of the window applied to the start —
        % sgolay's 'same' convolution is centred, so sample (win+1)/2 is
        % the first fully-supported derivative estimate.
        sIdx = (win+1)/2;
        v1 = dth1(sIdx);
        v2 = dth2(sIdx);
    end

    x0 = [th1(1); th2(1); v1; v2];
end


function plotFullTraceComparison(expData, fieldName, p_hat, FileName)
    % Simulate the full (unsegmented) trace with the fitted parameters
    % and overlay against measurement. This is the visual fit metric.

    for i = 1:numel(expData)
        d = expData(i).(fieldName);
        t  = d.t(:) - d.t(1);
        u  = d.u(:);
        th1_m = d.theta_1(:);
        th2_m = d.theta_2(:);

        % Initial state from data
        x0 = estimateInitialStateSG(t, th1_m, th2_m);

        % Flip input sign (matches iddata convention)
        u_flipped = -u;

        % Pack parameters into the call signature
        pa = num2cell(p_hat);
        u_fun = @(tt) zoh_input(tt, t, u_flipped);
        odefun = @(tt, xx) callRod1Model(FileName, tt, xx, u_fun(tt), pa);

        opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.01);
        [tSim, X] = ode45(odefun, t, x0, opts);

        figure('Name', sprintf('Full trace %s — %s', fieldName, expData(i).name), ...
               'NumberTitle', 'off');
        subplot(3,1,1);
        plot(t, u, 'k'); ylabel('u (raw)'); grid on;
        title(sprintf('Full %s trace: %s', fieldName, strrep(expData(i).name,'_','\_')));

        subplot(3,1,2); hold on; grid on;
        plot(t,    th1_m,  'Color',[.5 .5 .5], 'DisplayName','measured');
        plot(tSim, X(:,1), 'b',                'DisplayName','simulated');
        ylabel('\theta_1 [rad]'); legend('Location','best');

        subplot(3,1,3); hold on; grid on;
        plot(t,    th2_m,  'Color',[.5 .5 .5], 'DisplayName','measured');
        plot(tSim, X(:,2), 'b',                'DisplayName','simulated');
        ylabel('\theta_2 [rad]'); xlabel('t [s]'); legend('Location','best');
    end
end


function dx = callRod1Model(FileName, t, x, u, pa)
    fh = str2func(FileName);
    [dx, ~] = fh(t, x, u, pa{:});
end


function u_now = zoh_input(tt, t_grid, u_grid)
    if tt <= t_grid(1)
        u_now = u_grid(1);
    elseif tt >= t_grid(end)
        u_now = u_grid(end);
    else
        idx = find(t_grid <= tt, 1, 'last');
        u_now = u_grid(idx);
    end
end