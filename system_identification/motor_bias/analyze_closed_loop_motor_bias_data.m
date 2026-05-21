%% analyze_closed_loop_motor_bias_data.m
% Basic post-processing for the closed-loop motor bias/deadzone test.
%
% Edit the variable extraction section to match your logged signal names.
% Goal:
%   - Check whether positive and negative probe commands produce asymmetric response.
%   - Check whether the controller output required to hold angle has a systematic bias.

clear;
clc;
close all;

%% ================================================================
% Load your logged data
% ================================================================

% Example:
% load('your_logged_motor_bias_test.mat');
%
% Expected signals after you adapt this section:
%   t              [Nx1] time [s]
%   theta1_meas    [Nx1] measured theta_1 [rad]
%   theta1_ref     [Nx1] reference theta_1 [rad]
%   u_probe        [Nx1] probe command
%   u_cmd          [Nx1] actual motor command sent to motor

error('Edit this script first: load your logged data and define t, theta1_meas, theta1_ref, u_probe, and u_cmd.');

%% ================================================================
% Derivatives and residuals
% ================================================================

Ts = median(diff(t));
theta1_dot = gradient(theta1_meas, Ts);
theta1_acc = gradient(theta1_dot, Ts);

tracking_error = theta1_ref - theta1_meas;

%% ================================================================
% Plots
% ================================================================

figure('Name', 'Closed-loop motor bias test overview', ...
    'Units', 'normalized', 'Position', [0.06 0.08 0.88 0.82]);

tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t, rad2deg(theta1_meas), 'DisplayName', 'measured');
hold on;
plot(t, rad2deg(theta1_ref), '--', 'DisplayName', 'reference');
grid on;
ylabel('\theta_1 [deg]');
legend('Location', 'best');

nexttile;
plot(t, rad2deg(tracking_error));
grid on;
ylabel('error [deg]');

nexttile;
plot(t, u_probe, 'DisplayName', 'u probe');
hold on;
plot(t, u_cmd, 'DisplayName', 'u command');
grid on;
ylabel('input');
legend('Location', 'best');

nexttile;
plot(t, theta1_acc);
grid on;
ylabel('\alpha_1 [rad/s^2]');
xlabel('Time [s]');

figure('Name', 'Probe response asymmetry', ...
    'Units', 'normalized', 'Position', [0.16 0.16 0.7 0.55]);
scatter(u_probe, theta1_acc, 8, 'filled');
grid on;
xlabel('u_{probe}');
ylabel('estimated \alpha_1 [rad/s^2]');
title('Probe input versus estimated acceleration');

figure('Name', 'Actual motor command asymmetry', ...
    'Units', 'normalized', 'Position', [0.18 0.18 0.7 0.55]);
scatter(u_cmd, theta1_acc, 8, 'filled');
grid on;
xlabel('u_{cmd}');
ylabel('estimated \alpha_1 [rad/s^2]');
title('Actual motor command versus estimated acceleration');
