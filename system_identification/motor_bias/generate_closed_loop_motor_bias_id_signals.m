%% generate_closed_loop_motor_bias_id_signals.m
% Closed-loop motor bias / deadzone identification signals.
%
% Use these signals in Simulink:
%   theta1_ref_bias_id  -> reference input for a feedback controller
%   u_probe_bias_id     -> small probe command added inside the controller
%
% DO NOT feed u_probe_bias_id directly to the motor without feedback.
%
% Recommended Simulink structure:
%
%   From Workspace: theta1_ref_bias_id ----\
%                                      MATLAB Function controller ---> motor command
%   From Workspace: u_probe_bias_id --------/
%   measured theta_1 ----------------------/
%
% The MATLAB Function controller should implement saturation and angle safety.
% See motor_bias_feedback_controller.m.

clear;
clc;
close all;

%% ================================================================
% User settings
% ================================================================

Ts = 0.01;                         % sample time [s]

% Safe theta_1 range. The reference trajectory stays well inside this range.
thetaSoftLimitDeg = 70;            % controller enters recovery mode near this
thetaHardLimitDeg = 80;            % emergency limit / stop condition

% Reference schedule for theta_1 in degrees, relative to your chosen safe centre.
% Keep this conservative first. Increase only after checking the response.
refLevelsDeg = [0 10 -10 20 -20 30 -30 0];
refHoldTime = 8.0;                 % seconds per reference level
maxRefRateDegPerSec = 8.0;         % smooth reference slope limit

% Small bipolar probe amplitudes. These are added on top of the feedback command.
% They are intentionally small because feedback is already holding the system.
probeAmpList = [0.015 0.030 0.045 0.060];
probePulseTime = 0.6;              % seconds
probeRestTime  = 1.4;              % seconds

% Command safety for the generated probe signal only.
% The controller has its own final saturation.
maxAbsProbe = 0.075;

preRestTime  = 4.0;
postRestTime = 4.0;

%% ================================================================
% Build reference and probe sequence
% ================================================================

thetaRefDegRaw = [];
uProbeRaw = [];
segmentInfo = {};

% Initial rest at zero reference, no probe.
[thetaRefDegRaw, uProbeRaw, segmentInfo] = appendSegment( ...
    thetaRefDegRaw, uProbeRaw, segmentInfo, 0, 0, preRestTime, Ts, "initial rest");

for k = 1:numel(refLevelsDeg)
    refDeg = refLevelsDeg(k);

    % Build a probe sequence that fits inside the hold interval.
    [probeSegment, probeLabels] = makeProbeSegment( ...
        probeAmpList, probePulseTime, probeRestTime, refHoldTime, Ts, maxAbsProbe);

    n = numel(probeSegment);
    thetaSegment = refDeg * ones(n, 1);

    startIndex = numel(thetaRefDegRaw) + 1;
    thetaRefDegRaw = [thetaRefDegRaw; thetaSegment]; %#ok<AGROW>
    uProbeRaw      = [uProbeRaw; probeSegment]; %#ok<AGROW>
    endIndex = numel(thetaRefDegRaw);

    info.label = sprintf('reference %+g deg with probe sequence', refDeg);
    info.startIndex = startIndex;
    info.endIndex = endIndex;
    info.startTime = (startIndex - 1) * Ts;
    info.endTime = (endIndex - 1) * Ts;
    info.probeLabels = probeLabels;
    segmentInfo{end+1, 1} = info; %#ok<SAGROW>
end

% Final rest.
[thetaRefDegRaw, uProbeRaw, segmentInfo] = appendSegment( ...
    thetaRefDegRaw, uProbeRaw, segmentInfo, 0, 0, postRestTime, Ts, "final rest");

% Smooth/rate-limit reference to avoid aggressive position steps.
thetaRefDeg = rateLimitSignal(thetaRefDegRaw, Ts, maxRefRateDegPerSec);

% Convert to radians for Simulink controller.
thetaRefRad = deg2rad(thetaRefDeg(:));
uProbe = uProbeRaw(:);

t = (0:numel(thetaRefRad)-1)' * Ts;

%% ================================================================
% Create timeseries for From Workspace blocks
% ================================================================

theta1_ref_bias_id = timeseries(thetaRefRad, t);
theta1_ref_bias_id.Name = 'theta1_ref_bias_id';

u_probe_bias_id = timeseries(uProbe, t);
u_probe_bias_id.Name = 'u_probe_bias_id';

theta_soft_limit_rad = deg2rad(thetaSoftLimitDeg);
theta_hard_limit_rad = deg2rad(thetaHardLimitDeg);

assignin('base', 'theta1_ref_bias_id', theta1_ref_bias_id);
assignin('base', 'u_probe_bias_id', u_probe_bias_id);
assignin('base', 'theta_soft_limit_rad', theta_soft_limit_rad);
assignin('base', 'theta_hard_limit_rad', theta_hard_limit_rad);
assignin('base', 'segmentInfo', segmentInfo);

save('closed_loop_motor_bias_id_signals.mat', ...
    'theta1_ref_bias_id', ...
    'u_probe_bias_id', ...
    'theta_soft_limit_rad', ...
    'theta_hard_limit_rad', ...
    'segmentInfo', ...
    'Ts');

fprintf('\nGenerated closed-loop motor bias identification signals.\n');
fprintf('Total duration: %.2f s\n', t(end));
fprintf('Sample time: %.4f s\n', Ts);
fprintf('Reference range: %.1f to %.1f deg\n', min(thetaRefDeg), max(thetaRefDeg));
fprintf('Probe range: %.4f to %.4f\n', min(uProbe), max(uProbe));
fprintf('Saved to: closed_loop_motor_bias_id_signals.mat\n');
fprintf('\nUse From Workspace variables:\n');
fprintf('  theta1_ref_bias_id\n');
fprintf('  u_probe_bias_id\n');

%% ================================================================
% Plot
% ================================================================

figure('Name', 'Closed-loop motor bias identification schedule', ...
    'Units', 'normalized', 'Position', [0.08 0.12 0.84 0.72]);

tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t, rad2deg(thetaRefRad), 'LineWidth', 1.2);
grid on;
ylabel('\theta_{1,ref} [deg]');
title('Feedback reference trajectory');
yline(thetaSoftLimitDeg, '--');
yline(-thetaSoftLimitDeg, '--');
yline(thetaHardLimitDeg, ':');
yline(-thetaHardLimitDeg, ':');

nexttile;
plot(t, uProbe, 'LineWidth', 1.2);
grid on;
ylabel('u_{probe}');
xlabel('Time [s]');
title('Small probe command added inside feedback controller');
yline(0, '--');

%% ================================================================
% Local helper functions
% ================================================================

function [thetaRef, uProbe, segmentInfo] = appendSegment(thetaRef, uProbe, segmentInfo, refDeg, probeValue, duration, Ts, label)
    n = max(1, round(duration / Ts));
    startIndex = numel(thetaRef) + 1;
    thetaRef = [thetaRef; refDeg * ones(n, 1)]; %#ok<AGROW>
    uProbe   = [uProbe; probeValue * ones(n, 1)]; %#ok<AGROW>
    endIndex = numel(thetaRef);

    info.label = char(label);
    info.startIndex = startIndex;
    info.endIndex = endIndex;
    info.startTime = (startIndex - 1) * Ts;
    info.endTime = (endIndex - 1) * Ts;
    segmentInfo{end+1, 1} = info;
end

function [probeSegment, labels] = makeProbeSegment(ampList, pulseTime, restTime, totalTime, Ts, maxAbsProbe)
    probeSegment = [];
    labels = {};

    for k = 1:numel(ampList)
        A = min(abs(ampList(k)), maxAbsProbe);

        [probeSegment, labels] = appendProbeValue(probeSegment, labels, +A, pulseTime, Ts, sprintf('+%.3f pulse', A));
        [probeSegment, labels] = appendProbeValue(probeSegment, labels,  0, restTime,  Ts, 'rest');
        [probeSegment, labels] = appendProbeValue(probeSegment, labels, -A, pulseTime, Ts, sprintf('-%.3f pulse', A));
        [probeSegment, labels] = appendProbeValue(probeSegment, labels,  0, restTime,  Ts, 'rest');

        % Repeat in opposite order to reduce ordering effects.
        [probeSegment, labels] = appendProbeValue(probeSegment, labels, -A, pulseTime, Ts, sprintf('-%.3f pulse', A));
        [probeSegment, labels] = appendProbeValue(probeSegment, labels,  0, restTime,  Ts, 'rest');
        [probeSegment, labels] = appendProbeValue(probeSegment, labels, +A, pulseTime, Ts, sprintf('+%.3f pulse', A));
        [probeSegment, labels] = appendProbeValue(probeSegment, labels,  0, restTime,  Ts, 'rest');
    end

    nTotal = max(1, round(totalTime / Ts));

    if numel(probeSegment) < nTotal
        probeSegment = [probeSegment; zeros(nTotal - numel(probeSegment), 1)];
    else
        probeSegment = probeSegment(1:nTotal);
    end
end

function [x, labels] = appendProbeValue(x, labels, value, duration, Ts, label)
    n = max(1, round(duration / Ts));
    x = [x; value * ones(n, 1)]; %#ok<AGROW>
    labels{end+1, 1} = char(label);
end

function y = rateLimitSignal(x, Ts, maxRate)
    x = x(:);
    y = zeros(size(x));
    y(1) = x(1);
    maxStep = maxRate * Ts;

    for k = 2:numel(x)
        delta = x(k) - y(k-1);
        delta = max(min(delta, maxStep), -maxStep);
        y(k) = y(k-1) + delta;
    end
end
