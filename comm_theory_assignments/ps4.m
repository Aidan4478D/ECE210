% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE300 Communication Theory - Problem Set 3
% Copyright (C) 2025 Aidan Cusa <aidancusa@gmail.com>

clc; clear; close all; 

%% Question 2
amplitudes = [2, 10, 2];
bandwidths = [10, 2, 20];
kf = 2;

n = length(amplitudes);
delta_f = zeros(n,1);
mod_index = zeros(n,1);
carson_bw = zeros(n,1);
universal_bw = zeros(n,1);
band_type = strings(n,1);

for i = 1:n
    [delta_f(i), mod_index(i), band_type(i), carson_bw(i)] = get_fm_results(kf, amplitudes(i), bandwidths(i));
end

% inserting universal curve values based on inspection from chart in the
% notes => gets normalized bandwidth B_T / delta_f at a particular beta.
universal_bw(1) = 10 * delta_f(1); % norm bandwidth @ B = 0.4 approx. 10
universal_bw(2) = 3 * delta_f(2); % norm bandwidth @ B = 10 approx. 3
universal_bw(3) = 16 * delta_f(3); % norm bandwidth @ B = 0.2 approx. 16

Summary = table( (1:n).', amplitudes.', delta_f, mod_index, band_type, carson_bw, universal_bw, ...
                 'VariableNames', {'Case','amplitude','delta_f','mod. index','band type','CarsonBW (MHz)','UniversalBW (MHz)'});

disp(Summary);


%% Question 3
fc = 10; % carrier frequency
Ac = 2; % carrier amplitude
fm = 2; % message frequency
ka = 2; % amplitude sensitivity
fs = 100; % sampling frequency
t = 0:1/fs:19.99; % time vector for 2000 samples

Am = [0.4, 0.45, 0.6];

for i = 1:length(Am)
    m_t = Am(i) * cos(2 * pi * fm * t);
    s_t = Ac * (1 + ka*m_t) .* cos(2 * pi * fc * t);

    mod_index = max(ka * m_t) * 100;

    is_overmod = "not overmodulated";
    if mod_index > 100
        is_overmod = "overmodulated";
    end
    fprintf("Modulation index is: %d%%, thus the signal is %s\n", mod_index, is_overmod)

    subplot(3, 1, i);
    plot(t, s_t);
    title("AM Signal - A_m = " + Am(i))
    xlabel("Time (s)")
    ylabel("Amplitude (V)");
    grid on;
end


function [delta_f, mod_index, band_type, carson_bw] = get_fm_results(kf, amplitude, bandwidth)
    delta_f = kf * amplitude;
    mod_index = delta_f / bandwidth;
    if mod_index < 1
        band_type = 'Narrowband';
    else
        band_type = "Wideband";
    end
    carson_bw = 2 * (delta_f + bandwidth);
end
