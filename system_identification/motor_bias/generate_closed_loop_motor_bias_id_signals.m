%% generate_closed_loop_motor_bias_id_signals.m
% Generates closed-loop motor bias / deadzone identification signals.
%
% Use in Simulink:
%   From Workspace block 1: u_probe_bias_id
%   From Workspace block 2: theta1_ref_bias_id
%
% The actual motor command should be computed in a MATLAB Function block
% using theta_1 feedback and u_probe_bias_id.
%
% theta_1 equilibrium is assumed to be pi radians.
% Safety range is defined around theta_1 = pi:
%
%   theta_1 must stay close to pi +/- 80 degrees.

clear;
clc;
close all;

%% ================================================================
% User settings
% ================================================================

Ts = 0.01;                         % sample time [s]

% Stable hanging equilibrium of theta_1
theta1_equilibrium = pi;           % [rad]

% Safety limits around theta_1 = pi
thetaSoftLimitDeg = 70;            % soft recovery starts at +/-70 deg from pi
thetaHardLimitDeg = 80;            % hard limit at +/-80 deg from pi

theta_soft_limit_rad = deg2rad(thetaSoftLimitDeg);
theta_hard_limit_rad = deg2rad(thetaHardLimitDeg);

theta1_min_soft = theta1_equilibrium - theta_soft_limit_rad;
theta1_max_soft = theta1_equilibrium + theta_soft_limit_rad;

theta1_min_hard = theta1_equilibrium - theta_hard_limit_rad;
theta1_max_hard = theta1_equilibrium + theta_hard_limit_rad;

% Probe amplitudes.
% Start conservative. Increase only after checking the recorded theta_1.
ampList = [0.2 0.3 0.4 0.5 0.6 0.8 0.9];

% Timing settings
preRestTime  = 3.0;                % initial rest [s]
pulseTime    = 1.0;                % duration of each probe pulse [s]
restTime     = 2.0;                % rest between probe pulses [s]
postRestTime = 3.0;                % final rest [s]

% Number of repeats for each amplitude
nRepeats = 2;

% Input sign convention for the probe.
% Keep +1 unless you know your command sign should be reversed.
inputPolarity = 1;

% Maximum absolute probe command.
% The feedback controller will also saturate the final motor command.
maxAbsProbe = 1;

% Smooth command transitions slightly
useRampEdges = true;
rampTime = 0.15;                   % [s]

%% ================================================================
% Build probe input sequence
% ================================================================

u_probe = [];
segmentInfo = {};

[u_probe, segmentInfo] = appendConstant( ...
    u_probe, segmentInfo, 0, preRestTime, Ts, "initial rest");

for k = 1:numel(ampList)
    A = ampList(k);

    if A > maxAbsProbe
        warning('Amplitude %.4f exceeds maxAbsProbe %.4f. Clipping.', A, maxAbsProbe);
        A = maxAbsProbe;
    end

    for r = 1:nRepeats

        % Positive first
        [u_probe, segmentInfo] = appendPulsePair( ...
            u_probe, segmentInfo, +A, pulseTime, restTime, Ts, ...
            sprintf("A=%.3f repeat %d positive-first", A, r));

        % Negative first
        [u_probe, segmentInfo] = appendPulsePair( ...
            u_probe, segmentInfo, -A, pulseTime, restTime, Ts, ...
            sprintf("A=%.3f repeat %d negative-first", A, r));
    end
end

[u_probe, segmentInfo] = appendConstant( ...
    u_probe, segmentInfo, 0, postRestTime, Ts, "final rest");

% Apply sign convention and probe saturation
u_probe = inputPolarity * u_probe;
u_probe = max(min(u_probe, maxAbsProbe), -maxAbsProbe);

% Optional ramp smoothing
if useRampEdges
    u_probe = smoothCommandEdges(u_probe, Ts, rampTime);
end

% Time vector
t = (0:numel(u_probe)-1)' * Ts;
u_probe = u_probe(:);

%% ================================================================
% Create reference signal around theta_1 = pi
% ================================================================

theta1_ref = theta1_equilibrium * ones(size(t));

u_probe_bias_id = timeseries(u_probe, t);
u_probe_bias_id.Name = 'u_probe_bias_id';

theta1_ref_bias_id = timeseries(theta1_ref, t);
theta1_ref_bias_id.Name = 'theta1_ref_bias_id';

%% ================================================================
% Export to base workspace
% ================================================================

assignin('base', 'u_probe_bias_id', u_probe_bias_id);
assignin('base', 'theta1_ref_bias_id', theta1_ref_bias_id);

assignin('base', 'theta1_equilibrium', theta1_equilibrium);
assignin('base', 'theta_soft_limit_rad', theta_soft_limit_rad);
assignin('base', 'theta_hard_limit_rad', theta_hard_limit_rad);

assignin('base', 'theta1_min_soft', theta1_min_soft);
assignin('base', 'theta1_max_soft', theta1_max_soft);
assignin('base', 'theta1_min_hard', theta1_min_hard);
assignin('base', 'theta1_max_hard', theta1_max_hard);

assignin('base', 'segmentInfo', segmentInfo);
assignin('base', 'Ts_bias_id', Ts);

save('closed_loop_motor_bias_id_signals.mat', ...
    'u_probe_bias_id', ...
    'theta1_ref_bias_id', ...
    'theta1_equilibrium', ...
    'theta_soft_limit_rad', ...
    'theta_hard_limit_rad', ...
    'theta1_min_soft', ...
    'theta1_max_soft', ...
    'theta1_min_hard', ...
    'theta1_max_hard', ...
    'segmentInfo', ...
    'Ts');

fprintf('\nGenerated closed-loop motor bias identification signals.\n');
fprintf('Total duration: %.2f s\n', t(end));
fprintf('Sample time: %.4f s\n', Ts);
fprintf('Maximum absolute probe command: %.4f\n', max(abs(u_probe)));
fprintf('\nUse these From Workspace variables:\n');
fprintf('  u_probe_bias_id\n');
fprintf('  theta1_ref_bias_id\n');
fprintf('\nSafety range around theta_1 = pi:\n');
fprintf('  Soft range: %.4f rad to %.4f rad\n', theta1_min_soft, theta1_max_soft);
fprintf('  Hard range: %.4f rad to %.4f rad\n', theta1_min_hard, theta1_max_hard);
fprintf('\nSaved to closed_loop_motor_bias_id_signals.mat\n');

%% ================================================================
% Plot generated signals
% ================================================================

figure('Name', 'Closed-loop motor bias identification signals', ...
    'Units', 'normalized', 'Position', [0.1 0.15 0.8 0.65]);

tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t, u_probe, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('u_{probe}');
title('Probe input signal');

yline(0, '--');
ylim(1.2 * [-maxAbsProbe, maxAbsProbe]);

nexttile;
plot(t, theta1_ref, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('\theta_{1,ref} [rad]');
title('\theta_1 reference signal');

yline(theta1_equilibrium, '--', '\pi equilibrium');
yline(theta1_min_soft, '--', 'soft min');
yline(theta1_max_soft, '--', 'soft max');
yline(theta1_min_hard, ':', 'hard min');
yline(theta1_max_hard, ':', 'hard max');

%% ================================================================
% Local helper functions
% ================================================================

function [u, segmentInfo] = appendPulsePair(u, segmentInfo, A, pulseTime, restTime, Ts, label)

    [u, segmentInfo] = appendConstant( ...
        u, segmentInfo, A, pulseTime, Ts, label + " pulse 1");

    [u, segmentInfo] = appendConstant( ...
        u, segmentInfo, 0, restTime, Ts, label + " rest 1");

    [u, segmentInfo] = appendConstant( ...
        u, segmentInfo, -A, pulseTime, Ts, label + " pulse 2");

    [u, segmentInfo] = appendConstant( ...
        u, segmentInfo, 0, restTime, Ts, label + " rest 2");

end

function [u, segmentInfo] = appendConstant(u, segmentInfo, value, duration, Ts, label)

    n = max(1, round(duration / Ts));

    startIndex = numel(u) + 1;
    u = [u; value * ones(n, 1)];
    endIndex = numel(u);

    info.label = char(label);
    info.value = value;
    info.startIndex = startIndex;
    info.endIndex = endIndex;
    info.startTime = (startIndex - 1) * Ts;
    info.endTime = (endIndex - 1) * Ts;

    segmentInfo{end+1, 1} = info;

end

function uSmooth = smoothCommandEdges(u, Ts, rampTime)

    u = u(:);
    uSmooth = u;

    nRamp = max(1, round(rampTime / Ts));

    changeIdx = find(abs(diff(u)) > 1e-12) + 1;

    for k = 1:numel(changeIdx)
        idx0 = changeIdx(k);

        idxStart = idx0;
        idxEnd = min(numel(u), idx0 + nRamp - 1);

        if idx0 <= 1
            previousValue = u(idx0);
        else
            previousValue = uSmooth(idx0 - 1);
        end

        newValue = u(idx0);
        ramp = linspace(previousValue, newValue, idxEnd - idxStart + 1)';

        uSmooth(idxStart:idxEnd) = ramp;
    end

end