function [u_cmd, safety_mode, theta1_dot_est] = motor_bias_feedback_controller(theta1_meas, theta1_ref, u_probe, reset)
%MOTOR_BIAS_FEEDBACK_CONTROLLER Conservative feedback controller for motor bias tests.
%
% Use this code in a Simulink MATLAB Function block.
%
% Inputs:
%   theta1_meas : measured theta_1 [rad]
%   theta1_ref  : reference theta_1 [rad], from theta1_ref_bias_id
%   u_probe     : small probe command, from u_probe_bias_id
%   reset       : boolean reset signal; use false/0 if not needed
%
% Outputs:
%   u_cmd          : motor command
%   safety_mode    : 0 normal, 1 soft-limit recovery, 2 hard-limit recovery
%   theta1_dot_est : filtered velocity estimate [rad/s]
%
% Tune the constants below on the real setup.

%#codegen

%% ================================================================
% Tunable constants
% ================================================================

Ts = 0.01;                 % controller sample time [s]

% If positive controller action moves theta_1 the wrong way, set motorSign = -1.
motorSign = 1;

% Conservative PD gains. Start low and increase only if the system feels sluggish.
Kp = 0.45;
Kd = 0.08;

% Stronger gains used near the soft/hard angle limits.
KpSafety = 0.80;
KdSafety = 0.12;

% Command limits. Keep conservative for the first tests.
uMaxNormal = 0.18;
uMaxSafety = 0.25;

% Angle safety limits.
thetaSoftLimit = 70*pi/180;
thetaHardLimit = 80*pi/180;

% In safety mode, drive back toward this safe centre.
thetaRecoveryRef = 0.0;

% Low-pass filter for derivative estimate.
% Smaller alpha = smoother but more delayed.
alphaVel = 0.15;

%% ================================================================
% Persistent derivative estimate
% ================================================================

persistent thetaPrev thetaDotFilt initialized

if isempty(initialized) || reset ~= 0
    thetaPrev = theta1_meas;
    thetaDotFilt = 0.0;
    initialized = true;
end

rawDot = (theta1_meas - thetaPrev) / Ts;
thetaDotFilt = (1 - alphaVel) * thetaDotFilt + alphaVel * rawDot;
thetaPrev = theta1_meas;

theta1_dot_est = thetaDotFilt;

%% ================================================================
% Normal feedback law
% ================================================================

e = theta1_ref - theta1_meas;
u_fb = Kp * e - Kd * thetaDotFilt;

u_unsat = u_fb + u_probe;
u_limit = uMaxNormal;
safety_mode = 0;

%% ================================================================
% Safety override
% ================================================================

absTheta = abs(theta1_meas);

if absTheta >= thetaSoftLimit
    % Ignore probe and drive back to centre.
    eSafety = thetaRecoveryRef - theta1_meas;
    u_unsat = KpSafety * eSafety - KdSafety * thetaDotFilt;
    u_limit = uMaxSafety;
    safety_mode = 1;
end

if absTheta >= thetaHardLimit
    % Strong recovery. You can also connect safety_mode == 2 to a Stop block.
    eSafety = thetaRecoveryRef - theta1_meas;
    u_unsat = 1.25 * KpSafety * eSafety - 1.25 * KdSafety * thetaDotFilt;
    u_limit = uMaxSafety;
    safety_mode = 2;
end

%% ================================================================
% Final saturation and sign convention
% ================================================================

u_cmd = motorSign * saturate(u_unsat, -u_limit, u_limit);

end

function y = saturate(x, xmin, xmax)
    y = min(max(x, xmin), xmax);
end
