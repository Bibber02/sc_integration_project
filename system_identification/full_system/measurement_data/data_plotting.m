clear;
clc;
close all;

%% ================================================================
% Plot cleaned full-system identification data
%
% Each file should contain:
%   theta_1
%   theta_2
%   u_ts
%
% Folder structure:
%   prbs/
%   chirp/
% ================================================================

%% Settings

scriptFolder = fileparts(mfilename('fullpath'));
dataFolder = scriptFolder;

amplitudes = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];

%% Plot PRBS measurements

figure('Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);
tiledlayout(10, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:10

    ampText = strrep(sprintf('%.2f', amplitudes(k)), '.', 'p');
    filename = fullfile(dataFolder, 'prbs', sprintf('fullsystem_prbs_A%s_run%02d.mat', ampText, k));

    load(filename, 'theta_1', 'theta_2', 'u_ts');

    nexttile;
    plot(theta_1.Time, squeeze(theta_1.Data), 'LineWidth', 1.1);
    grid on;
    ylabel(sprintf('Run %02d', k));
    title('\theta_1');
    xlabel('Time [s]');

    nexttile;
    plot(theta_2.Time, squeeze(theta_2.Data), 'LineWidth', 1.1);
    grid on;
    title('\theta_2');
    xlabel('Time [s]');

    nexttile;
    plot(u_ts.Time, squeeze(u_ts.Data), 'LineWidth', 1.1);
    grid on;
    title('Input u');
    xlabel('Time [s]');
    ylim([-0.4 0.4]);

end

sgtitle('PRBS identification measurements');


%% Plot chirp measurements

figure('Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);
tiledlayout(10, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:10

    ampText = strrep(sprintf('%.2f', amplitudes(k)), '.', 'p');
    filename = fullfile(dataFolder, 'chirp', sprintf('fullsystem_chirp_A%s_run%02d.mat', ampText, k));

    load(filename, 'theta_1', 'theta_2', 'u_ts');

    nexttile;
    plot(theta_1.Time, squeeze(theta_1.Data), 'LineWidth', 1.1);
    grid on;
    ylabel(sprintf('Run %02d', k));
    title('\theta_1');
    xlabel('Time [s]');

    nexttile;
    plot(theta_2.Time, squeeze(theta_2.Data), 'LineWidth', 1.1);
    grid on;
    title('\theta_2');
    xlabel('Time [s]');

    nexttile;
    plot(u_ts.Time, squeeze(u_ts.Data), 'LineWidth', 1.1);
    grid on;
    title('Input u');
    xlabel('Time [s]');
    ylim([-0.4 0.4]);

end

sgtitle('Chirp identification measurements');
