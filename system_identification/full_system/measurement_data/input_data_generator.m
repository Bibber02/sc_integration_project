clear;
clc;
close all;

%% ================================================================
% Generate From Workspace input signals for full-system ID
%
% Each signal contains:
%   rest before excitation
%   20 s excitation
%   rest after excitation
%
% Generated variables:
%   u_prbs_01, ..., u_prbs_05
%   u_chirp_01, ..., u_chirp_05
%
% Also saved next to this script:
%   generated_id_inputs.mat
% ================================================================

%% Settings
scriptFolder = fileparts(mfilename('fullpath'));
generatedInputFile = fullfile(scriptFolder, 'generated_id_inputs.mat');

Ts = 0.01;                 % sample time [s], same as sensors
restBefore = 4.0;          % seconds
excitationTime = 20.0;     % seconds
restAfter = 4.0;           % seconds

% Adjust these if the motion is too small/large.
% If amplitudes below 0.2 do not move the motor, increase the lower values.
prbsAmplitudes  = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];
chirpAmplitudes = [0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34];

% PRBS switching time. Smaller = richer excitation, but more aggressive.
% Use 0.25 to 0.70 s as a reasonable range.
prbsSwitchTime = 0.25;     % seconds

% Chirp frequency range.
% Do not sweep too wide in 20 seconds.
chirpF0 = 0.2;             % Hz
chirpF1 = 3.0;             % Hz

% Optional final safety saturation.
uMax = 0.30;

rng(23);                    % repeatable PRBS signals

%% Time vector
tTotal = restBefore + excitationTime + restAfter;
t = (0:Ts:tTotal).';

idxExc = t >= restBefore & t < restBefore + excitationTime;
tExc = t(idxExc) - restBefore;

%% Generate PRBS signals
prbsNames = strings(length(prbsAmplitudes), 1);

for k = 1:length(prbsAmplitudes)

    A = prbsAmplitudes(k);

    u = zeros(size(t));

    % Generate block-wise random +/- 1 values during excitation.
    samplesPerSwitch = max(1, round(prbsSwitchTime / Ts));
    nExc = sum(idxExc);
    nBlocks = ceil(nExc / samplesPerSwitch);

    blockValues = 2 * randi([0 1], nBlocks, 1) - 1;   % +/- 1
    prbsExc = repelem(blockValues, samplesPerSwitch);
    prbsExc = prbsExc(1:nExc);

    u(idxExc) = A * prbsExc;

    % Safety saturation.
    u = max(-uMax, min(uMax, u));

    varName = sprintf('u_prbs_%02d', k);
    assignin('base', varName, timeseries(u, t));
    eval(sprintf('%s = timeseries(u, t);', varName));

    prbsNames(k) = varName;
end

%% Generate chirp signals
chirpNames = strings(length(chirpAmplitudes), 1);

for k = 1:length(chirpAmplitudes)

    A = chirpAmplitudes(k);

    u = zeros(size(t));

    % Linear chirp phase:
    % f(t) = f0 + (f1-f0)*t/T
    % phase = 2*pi*(f0*t + 0.5*(f1-f0)/T*t^2)
    T = excitationTime;
    phase = 2*pi*(chirpF0*tExc + 0.5*(chirpF1 - chirpF0)/T .* tExc.^2);

    chirpExc = A * sin(phase);

    u(idxExc) = chirpExc;

    % Safety saturation.
    u = max(-uMax, min(uMax, u));

    varName = sprintf('u_chirp_%02d', k);
    assignin('base', varName, timeseries(u, t));
    eval(sprintf('%s = timeseries(u, t);', varName));

    chirpNames(k) = varName;
end

%% Plot all signals
figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);
tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 8:10
    nexttile;
    sig = eval(prbsNames(k));
    plot(sig.Time, sig.Data);
    grid on;
    ylim([-uMax*1.1, uMax*1.1]);
    title(sprintf('%s, A = %.2f', prbsNames(k), prbsAmplitudes(k)));
    xlabel('Time [s]');
    ylabel('u');

    nexttile;
    sig = eval(chirpNames(k));
    plot(sig.Time, sig.Data);
    grid on;
    ylim([-uMax*1.1, uMax*1.1]);
    title(sprintf('%s, A = %.2f', chirpNames(k), chirpAmplitudes(k)));
    xlabel('Time [s]');
    ylabel('u');
end

sgtitle('Generated ID input signals');

%% Save all generated signals
metadataNames = {
    'Ts'
    'restBefore'
    'excitationTime'
    'restAfter'
    'prbsAmplitudes'
    'chirpAmplitudes'
    'prbsSwitchTime'
    'chirpF0'
    'chirpF1'
};

saveVars = [
    cellstr(prbsNames(:))
    cellstr(chirpNames(:))
    metadataNames
];

save(generatedInputFile, saveVars{:});

fprintf('\nGenerated variables:\n');
disp(prbsNames);
disp(chirpNames);

fprintf('\nSaved to %s\n', generatedInputFile);
