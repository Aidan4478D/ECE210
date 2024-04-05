% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment07.m -- Under Pressure
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 

fs = 80e3;   % sampling
fn = fs / 2; % nyquist

%% part 1
brake_data = readmatrix("40p_1000ms.csv");

%% part 2

t = brake_data(:, 1) ./ fs;
norm_brake_data = normalize(brake_data(:, 2), 'range', [0, 1]);

x = fft(norm_brake_data);
x = x(1:numel(x)/2); % "removing" negative frequencies by /2
f = fs * linspace(0, pi, numel(x)) / (2 * pi); % normalize fs to w range

fMaskLow = f >= 5; % only care about signals >= 5 hZ
x = x(fMaskLow);
f = f(fMaskLow);

dB = 20 * log10(abs(x));

fMaskHigh = f <= 40; % signals <= 40 hZ 
f = f(fMaskHigh);
dB = dB(fMaskHigh);

% plot raw dft
figure;
plot(f, dB)
title("Raw Brake Data");
yline(max(dB) - 20, '--', '-20dB');
yline(max(dB) - 40, '--', '-40dB');
xlabel('Frequency [hZ]');
ylabel('Magnitude [dB]');

%% part 3
Wp = 10 / fn; % ~10 Hz goes below max - 20dB
Ws = 20 / fn; % ~20 Hz goes below max - 40dB
Rp = 0.1; % passband ripple
Rs = 40;  % stopband ripple

% thus corner frequencies = 10 Hz, 20 Hz

[n, Wn] = ellipord(Wp, Ws, Rp, Rs);
[z, p, k] = ellip(n, Rp, Rs, Wp); 

%% part 4
[sos, gain] = zp2sos(z, p, k); % get sos digital filter
sos_range = linspace(0, Ws * fn + 20, 1e3); % plot up to fstop + 20
sos_response = gain * freqz(sos, sos_range, fs);
sos_response_dB = 20 * log10(abs(sos_response));

% plot filter
figure;
plot(sos_range, sos_response_dB); 
xline(Wp * fn, '--', 'Passband');
xline(Ws * fn, '--', 'Stopband');
yline(-Rp, '--', 'Passband');
yline(-Rs, '--', 'Stopband');
title("SOS Digital Filter");
xlabel('Frequency [hZ]');
ylabel('Magnitude [dB]');
ylim([-50, 10]);

%% part 5
S = gain * sosfilt(sos, norm_brake_data);

% plot raw data
figure;
plot(t, norm_brake_data);
title("Raw Data");
xlabel("Time [s]");

% plot filtered data
figure;
plot(t, S); 
ylim([0, 1]);
yline(0.33, '--', 'DC Bias'); %is this the correct DC bias?
title("Filtered Data");
xlabel("Time [s]");