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

    mod_index = max(ka * m_t) * 100; % value is in percent

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

%% Question 4

% set seed for reproducibility
rng(7)

% part C
sig2 = 1;
num_samples = 1e6;
threshold = 3.71692219; % computed by hand

nI = sqrt(sig2) * randn(num_samples, 1);
nQ = sqrt(sig2) * randn(num_samples, 1);
R = abs(nI + 1j*nQ);

plot_histogram(R, threshold, "C")

% part D
fraction_d = sum(R > threshold) / num_samples;
theoretical_prob_d = exp(-threshold^2 / (2*sig2));
error_d = abs(fraction_d - theoretical_prob_d)/theoretical_prob_d * 100;

fprintf("Part D:\n\tFraction of time R > threshold: %f\n", fraction_d);
fprintf("\tTheoretical probability: %f\n", theoretical_prob_d);
fprintf("\tRelative error %f%%\n", error_d);

% part E
sig2_min1dB = 10^(-0.1); % -1 dB
nI_e = sqrt(sig2_min1dB) * randn(num_samples, 1);
nQ_e = sqrt(sig2_min1dB) * randn(num_samples, 1);
R_e = abs(nI_e + 1j*nQ_e);

fraction_e = sum(R_e > threshold) / num_samples;
theoretical_prob_e = exp(-threshold^2 / (2*sig2_min1dB));
error_e = abs(fraction_e - theoretical_prob_e)/theoretical_prob_e * 100;

fprintf("Part E:\n\tFraction of time R > threshold: %f\n", fraction_e);
fprintf("\tTheoretical probability: %f\n", theoretical_prob_e);
fprintf("\tRelative error %f%%\n", error_e);

plot_histogram(R_e, threshold, "E")

% Part F
sig2_plus1dB = 10^(0.1); % +1 dB
nI_f = sqrt(sig2_plus1dB) * randn(num_samples, 1);
nQ_f = sqrt(sig2_plus1dB) * randn(num_samples, 1);
R_f = abs(nI_f + 1j*nQ_f);

fraction_f = sum(R_f > threshold) / num_samples;
theoretical_prob_f = exp(-threshold^2 / (2*sig2_plus1dB));
error_f = abs(fraction_f - theoretical_prob_f)/theoretical_prob_f * 100;

fprintf("Part F:\n\tFraction of time R > threshold: %f\n", fraction_f);
fprintf("\tTheoretical probability: %f\n", theoretical_prob_f);
fprintf("\tRelative error %f%%\n", error_f);

plot_histogram(R_f, threshold, "F")

% Part G:
    % Yes, the probability of exceeding the threshold seems particularly
    % sensitive to the power level. A seeminly small +-1 dB change in the power  
    % Comparing the magnitude of probability of exceeding the threshold between
    % the baseline (part D), -1dB (part E), and +1dB (part F), we find:
    
    % -1 dB compared to baseline => 1.590e-4 / 1.026e-3 = 0.155 (6.5x smaller)
    % +1 dB compared to baseline => 4.128e-3 / 1.026e-3 = 4.03 (4x larger) 
    % +1 dB compared to -1 dB => 4.128e-3 / 1.590e-4 = 25.96 (26x larger, 
    % over an order of magnitude)
    
    % This occurs becasue for the Rayleigh pdf, the tail probability 
    % depends exponentially on 1/sigma^2, so small dB changes in sigma^2 will 
    % cause large multiplicative changes in P(R > p).

% Part H:
threshold_H = sqrt(-2 * sig2 * log(1e-4));

theo_H_base = exp(-threshold_H^2 / (2*sig2));
theo_H_min1dB = exp(-threshold_H^2 / (2*sig2_min1dB));
theo_H_plus1dB = exp(-threshold_H^2 / (2*sig2_plus1dB));

fprintf('\nPart H (theoretical values for P(R>p)=10^-4):\n');
fprintf('\tBaseline (0 dB): P(R>p) = %f\n', theo_H_base);
fprintf('\t-1 dB: P(R>p) = %f\n', theo_H_min1dB);
fprintf('\t+1 dB: P(R>p) = %f\n', theo_H_plus1dB);

% calculate ratios
min1_over_base = theo_H_min1dB / theo_H_base;
plus1_over_base = theo_H_plus1dB / theo_H_base;
plus1_over_min1 = theo_H_plus1dB / theo_H_min1dB;

fprintf('\t-1 dB / base ratio: %f\n', min1_over_base);
fprintf('\t+1 dB / base ratio: %f\n', plus1_over_base);
fprintf('\t+1 dB / -1 dB ratio: %f\n', plus1_over_min1);

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

function plot_histogram(R, threshold, q_letter)
    figure;
    histogram(R, 100); 
    hold on;
    xline(threshold, 'r-', 'LineWidth', 1);
    title("Q4 - (Part " + q_letter + ")")
    xlabel('R');
    ylabel('Density');
    legend('Simulated R', sprintf('Threshold = %.2f', threshold));
    grid on;
end