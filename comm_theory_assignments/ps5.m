% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE300 Communication Theory - Problem Set 5
% Copyright (C) 2025 Aidan Cusa <aidancusa@gmail.com>

clc; clear; close all; 

%% Q1

snr = 4; % 4 dB / bit
gamma = 10^(snr / 10);

% orthogonal case => p=0
BER_orth = qfunc(sqrt(gamma));

% lowest BER case => p=-0.217
BER_low = qfunc(sqrt((1 - 0.217) * gamma));

fprintf("Question 1:\n\tOrthogonal case: %f\n", BER_orth);
fprintf("\tLowest BER case: %f\n", BER_low);

% part B
% non-coherent detection, equation from chart on deck 04 slide 17
BER_noncoh_orth = 0.5 * exp(-gamma/2);
fprintf('\tnoncoherent orthogonal BER: %.6e\n\n', BER_noncoh_orth);

% need to find BER for noncoherent detection for optimal case


%% Q2

% setting E0 = 1 for ease of calculations
constellation = [-3 - 1j, -3 + 1j, ...
                 -1 - 1j, -1 + 1j, ...
                  3 - 1j, 3 + 1j, ...
                  1 - 1j, 1 + 1j];

% part B
%display(constellation)

E_s = mean(abs(constellation).^2);
E_b = E_s / 3; % 3 bits / sym
fprintf("Question 2\nE_b = %f * E_0\n", E_b);

% part C
M = length(constellation);
distances = [];
for i = 1:M
    for j = 1:M
        if i ~= j
            d = abs(constellation(i) - constellation(j));
            distances = [distances; d];
        end
    end
end

unique_distances = unique(distances);
%display(unique_distances)

fprintf('part c: distance values\n');
for k = 1:length(unique_distances)
    d_normalized = unique_distances(k);
    count = sum(abs(distances - unique_distances(k)) < 1e-6); % within threshold to account for imprecise storage
    fprintf('\td_%d = %f * sqrt(E0), count = %d pairs\n', k, d_normalized, count);
end
fprintf('\ttotal pairs: %d (should 56)\n\n', sum(abs(distances) > 0));

% part D

% equation from deck 02 slide 20 
fprintf('part d: P_e <= sum(a_k * Q(sqrt(b_k * gamma_b)))\n')
a_k = [];
b_k = [];
for k = 1:length(unique_distances)
    d_k = unique_distances(k);
    count = sum(abs(distances - d_k) < 1e-6);

    a_k(k) = count / (M);
    b_k(k) = (d_k^2) / (2 * E_b);

    fprintf('\tk=%d: d_k=%f * sqrt(E0), a_k=%f, b_k=%f\n', k, d_k, a_k(k), b_k(k));
end

% part E
gamma_b_dB = 0:0.5:20;
gamma_b_lin = 10 .^ (gamma_b_dB / 10);

pe_detailed = zeros(size(gamma_b_lin));
pe_simplified = zeros(size(gamma_b_lin));

for idx = 1:length(gamma_b_lin)
    gb = gamma_b_lin(idx);

    % Detailed union bound
    for k = 1:length(a_k)
        pe_detailed(idx) = pe_detailed(idx) + a_k(k) * qfunc(sqrt(b_k(k) * gb));
    end

    % Simplified (dominant term with d_min only)
    pe_simplified(idx) = get_pe_estimate(gb);
end

% semilogy to plot with a logarithmic scale for the y-axis and a linear scale for the x-axis
figure('Position', [100, 100, 800, 600]);
semilogy(gamma_b_dB, pe_detailed, 'b-', 'LineWidth', 2);
hold on;
semilogy(gamma_b_dB, pe_simplified, 'r--', 'LineWidth', 2);
grid on;
xlabel('\gamma_b (dB/bit)');
ylabel('p_e');
title('8-Point Constellation - Symbol Error Probability');
legend('Detailed Union Bound', 'Simplified (d_{min} only)');
ylim([1e-6 1e-2]);


%% Q4

% part d
snr_db = 15;
snr_linear = 10^(snr_db/10);

signal_power = 1;
I_worst = 0.12; % 0.06 + 0.04 + 0.02

SIR_dB = 10 * log10 (signal_power / I_worst^2);

noise_power = signal_power / snr_linear;
SNIR_dB = 10 * log10 (signal_power / (I_worst^2 + noise_power));

fprintf("\nQuestion 4\npart d:\n")
fprintf("\tSIR (dB) = %f\n", SIR_dB);
fprintf("\tSNIR (dB) = %f\n", SNIR_dB);

%% Q6
fprintf("\nQuestion 6\n")
R_s = 48e3; % symbol rate
beta = 0.2; % rolloff
sps = 16; % samples per sym

% part a
bandwidth = (1 + beta) * R_s / 2;
fprintf("part a:\n\tbandwidth=%f\n", 2 * bandwidth); % 2x for two-sided bandwidth

% part b
fs_digital = sps * R_s;
fprintf("part b:\n\tdigital sampling rate=%f\n", fs_digital);

% part c
bit_rate = 2 * R_s; % QPSK = log2(4) = 2 bits/sym
fprintf("part c:\n\tbit rate=%f\n", bit_rate);

% part d
span = 3;

rrc = rcosdesign(beta, span, sps, 'sqrt');
mf = conv(rrc, fliplr(rrc));  % matched filter output

% plot pulse shapes
figure('Position', [100, 100, 1200, 400]);
subplot(2,1,1);
stem(rrc, 'filled', 'MarkerSize', 3);
grid on;
title('Transmit RRC Pulse Shape');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(2,1,2);
stem(mf, 'filled', 'MarkerSize', 3);
grid on;
title('Matched Filter Output');
xlabel('Sample Index');
ylabel('Amplitude');

% part e 

% get peak index
[~, center_idx] = max(abs(mf));

% get samples around peak
m_range = -span:span;
gd = zeros(size(m_range));
for idx = 1:length(m_range)
    m = m_range(idx);
    gd(idx) = mf(center_idx + m*sps);
end

%display(gd)

% normalize so gd[0] = 1
gd0 = gd(m_range==0);
gd = gd / gd0;

A_max = 1;  % QPSK symbols have unit magnitude
I = 0;
for i = 1:length(m_range)
    if m_range(i) ~= 0
        I = I + A_max * abs(gd(i));
    end
end

signal_power = (A_max)^2;
SIR0_dB = 10*log10(signal_power / I^2);

fprintf("part e:\n\tSIR0 = %.4f dB\n", SIR0_dB);

% part f
SIR0_lin = 10^(SIR0_dB/10);

SNIR_target_dB  = SIR0_dB - 5;
SNIR_target_lin = 10^(SNIR_target_dB/10);

% SNR = SNIR * SIR0 / (SIR0 - SNIR)
SNR_required_lin = (SNIR_target_lin * SIR0_lin) / (SIR0_lin - SNIR_target_lin);
SNR_required_dB  = 10*log10(SNR_required_lin);

fprintf("part f:\n\ttarget SNIR = %.4f dB\n", SNIR_target_dB);
fprintf("\trequired SNR = %.4f dB\n", SNR_required_dB);


% idk why but my MATLAB is not having it with "qfunc" so I just made my own
function [value] = qfunc(x)
    value = 0.5*erfc(x/sqrt(2));
end

% function that computes estimates of P_e as a function of gamma for d_min
% coefficients (a = 2.5, b = 1)
function [p_e] = get_pe_estimate(gamma_k)
    p_e = 2.5 * qfunc(sqrt(gamma_k));
end















