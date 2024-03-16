% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment06.m -- Four-Year Transform
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 

%% part 1
db2mag_ = @(x) 10 .^ (x / 20);

%% part 2
fs = 96e3;
n = 192e3;

freqs = [-20.48e3; -360; 996; 19.84e3];
H = db2mag_([14; -10; 0; 2]);

T = 1/fs;
t = (0:n - 1) * T;

s = sum(H .* exp(1j * 2 * pi * freqs .* t)); % superposition of h*e^(2pift)
noise = db2mag_(-10) * randn(size(s));

S = fftshift(fft(s + noise)); % generate signal

figure;
plot(fs / n * (-n / 2:n / 2 - 1), 20 * log10(abs(S)));
xlabel('Frequency [Hz]');
ylabel('Phase [Degrees]');
title("Phase Response");

%% part 3
k = 0.53; % scaling factor

zeros = [0.76 + 0.64j
         0.69 + 0.71j
         0.82 + 0.57j];

% i wont lie i manually typed out every one then saw your solution 
% (after finished i promise 🙏) and was like damn that's kinda a good idea

zeros = [zeros; conj(zeros)]; % get for +- zeros

poles = [0.57 + 0.78j
         0.85 + 0.48j];

poles = [poles; conj(poles); 0.24; 0.64];

% get transfer function given zpk
[b, a] = zp2tf(zeros, poles, k);

% only care about frequencies up to nyquist bandwidth due to aliasing
w = linspace(0, fs / 2, 1e4);

h = freqz(b, a, w, fs); 
mag = 20 * log10(abs(h));
phase = unwrap(angle(h)) * 180/pi; 

% plot mag
figure;
subplot(2, 2, 2);
plot(w, mag);

grid on;
xlim([0, fs/2]);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title("Magnitude Response");

% plot phase
subplot(2, 2, 4);
plot(w, phase);

grid on;
xlim([0, fs/2]);
xlabel('Frequency [Hz]');
ylabel('Phase [Degrees]');
title("Phase Response");

% plot z-p
subplot(2, 2, [1,3]);
zplane(zeros, poles); 

sgtitle("Response of Filter H");

