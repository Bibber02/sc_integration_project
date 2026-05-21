function [u_cmd, safety_mode] = motor_bias_simple_feedback(theta1, u_probe)
% Simple closed-loop safety controller for motor bias/deadzone testing.
%
% The controller keeps theta_1 close to theta_1 = pi.
% Safety limits are based on the deviation from pi.
%
% Inputs:
%   theta1  = measured rod-1 angle [rad]
%   u_probe = small probe signal from From Workspace
%
% Outputs:
%   u_cmd       = actual command sent to motor
%   safety_mode = 0 normal, 1 soft recovery, 2 hard recovery

% ---------- tunable constants ----------
Kp = 0.45;

theta_ref = pi;              % stable hanging equilibrium [rad]

maxCommand = 0.15;           % saturation of actual motor command

softLimit = 70*pi/180;       % soft recovery starts at +/-70 deg from pi
hardLimit = 80*pi/180;       % hard recovery starts at +/-80 deg from pi

motorSign = 1;               % change to -1 if feedback pushes the wrong way

% ---------- angle error ----------
% Use this if theta_1 is continuous around pi:
theta_error = theta_ref - theta1;

% Use this instead if theta_1 wraps between -pi and pi:
% theta_error = atan2(sin(theta_ref - theta1), cos(theta_ref - theta1));

absThetaError = abs(theta_error);

% ---------- safety logic ----------
if absThetaError >= hardLimit

    safety_mode = 2;

    % Hard recovery: ignore probe and drive toward theta_1 = pi
    u_feedback = Kp * theta_error;
    u_cmd_raw = motorSign * u_feedback;

elseif absThetaError >= softLimit

    safety_mode = 1;

    % Soft recovery: ignore probe and drive toward theta_1 = pi
    u_feedback = Kp * theta_error;
    u_cmd_raw = motorSign * u_feedback;

else

    safety_mode = 0;

    % Normal test mode: feedback plus small probe signal
    u_feedback = Kp * theta_error;
    u_cmd_raw = motorSign * (u_feedback + u_probe);

end

% ---------- saturation ----------
if u_cmd_raw > maxCommand
    u_cmd = maxCommand;
elseif u_cmd_raw < -maxCommand
    u_cmd = -maxCommand;
else
    u_cmd = u_cmd_raw;
end