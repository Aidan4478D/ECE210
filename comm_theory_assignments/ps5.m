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
BER_low = qfunc(sqrt((1 - (-0.217)) * gamma));

fprintf("Question 1:\n\tOrthogonal case: %f\n", BER_orth);
fprintf("\tLowest BER case: %f\n", BER_low);

% part B
% values based on inspection from the chart in the slides
fprintf('\tnoncoherent orthogonal BER: 0.2\n');
fprintf('\tnoncoherent lowest BER: 0.02\n');

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

    % detailed union bound
    for k = 1:length(a_k)
        pe_detailed(idx) = pe_detailed(idx) + a_k(k) * qfunc(sqrt(b_k(k) * gb));
    end

    % simplified (dominant term with d_min only)
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
xlim([1 20]);


%% Q4

% part d
snr_db = 15;
snr_linear = 10^(snr_db/10);

signal_power = 1;
I_worst = 0.12; % 0.06 + 0.04 + 0.02

SIR_dB = 10 * log10 (signal_power / I_worst^2);

noise_power = signal_power / snr_linear;
% display(noise_power)

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
bandwidth = (1 + beta) * R_s / 2; % (one sided bandwidth)
fprintf("part a:\n\tbandwidth=%f\n", bandwidth);

% part b
fs_digital = sps * R_s;
fprintf("part b:\n\tdigital sampling rate=%f\n", fs_digital);

% part c
bit_rate = 2 * R_s; % QPSK = log2(4) = 2 bits/sym
fprintf("part c:\n\tbit rate=%f\n", bit_rate);

% part d
span = 3;

rrc = rcosdesign(beta, span, sps, 'sqrt');
mf = rrc; % mf = transmit filter for rrc

rrc_conv = conv(rrc, mf);  % matched filter output

% plot pulse shapes
figure('Position', [100, 100, 1200, 400]);
subplot(2,1,1);
stem(rrc, 'filled', 'MarkerSize', 3);
grid on;
title('Transmit RRC Pulse Shape');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(2,1,2);
stem(rrc_conv, 'filled', 'MarkerSize', 3);
grid on;
title('Matched Filter Output');
xlabel('Sample Index');
ylabel('Amplitude');

% part e
% get peak index
[~, center_idx] = max(abs(rrc_conv));

% get samples around peak
m_range = -span:span;
gd = zeros(size(m_range));
for idx = 1:length(m_range)
    m = m_range(idx);
    gd(idx) = rrc_conv(center_idx + m*sps);
end

%display(gd)

% find the center sample
gd0 = gd(m_range == 0);

A_max = 1;
signal_power = abs(gd0)^2;

I = A_max * sum(abs(gd(m_range ~= 0)));
SIR0_linear = signal_power / (I^2);
SIR0_dB = 10*log10(SIR0_linear);
fprintf("part e:\n\tSIR0 = %.4f dB\n", SIR0_dB);

% part f
SIR0_lin = 10^(SIR0_dB/10);

SNIR_target_dB  = SIR0_dB - 5;
SNIR_target_lin = 10^(SNIR_target_dB/10);

% SNR = SNIR * SIR0 / (SIR0 - SNIR)
SNR_required_lin = (SNIR_target_lin * SIR0_lin) / (SIR0_lin - SNIR_target_lin);
SNR_required_dB = 10*log10(SNR_required_lin);

fprintf("part f:\n\ttarget SNIR = %.4f dB\n", SNIR_target_dB);
fprintf("\trequired SNR = %.4f dB\n", SNR_required_dB);


%% Q7

fprintf("\nQuestion 7\n");

num_bits = 1e5;
bits = randi([0, 1], num_bits, 1);

% % QPSK constellation with gray coding (normalized to unit magnitude)
% 00 -> ( +1 + j)/sqrt(2)
% 01 -> ( -1 + j)/sqrt(2)
% 11 -> ( -1 - j)/sqrt(2)
% 10 -> ( +1 - j)/sqrt(2)

bit_pairs = reshape(bits, 2, []).';
sym = zeros(size(bit_pairs,1),1);

for n = 1:size(bit_pairs,1)
    b1 = bit_pairs(n,1); 
    b2 = bit_pairs(n,2);

    if b1==0 && b2==0
        sym(n) = (1 + 1j)/sqrt(2);
    elseif b1==0 && b2==1
        sym(n) = (-1 + 1j)/sqrt(2);
    elseif b1==1 && b2==1
        sym(n) = (-1 - 1j)/sqrt(2);
    elseif b1==1 && b2==0
        sym(n) = (1 - 1j)/sqrt(2);
    end
end

% upsample and pulse shape
sym_up = upsample(sym, sps);
x = filter(rrc, 1, sym_up);

% plot zoomed envelope so it doesn't look solid
num_syms_to_plot = 300;
Nplot = num_syms_to_plot * sps;

figure;
plot(abs(x(1:Nplot)), 'LineWidth', 1);
grid on;
title('Q7(a) Transmitted waveform envelope |x[n]| (first 300 symbols)');
xlabel('Sample index');
ylabel('Envelope');

% find average SIR in transmitted waveform
Lh = length(rrc);
d_tx = span*sps / 2; % group delay

idxs = d_tx+1 : sps : length(x);   % valid peak locations in x
Ns = min(length(sym), length(idxs));  % how many we can actually sample

x_samp = x(idxs(1:Ns));
sym_use = sym(1:Ns);

% ISI only
err_tx = x_samp - sym_use;

S_avg = mean(abs(sym_use).^2);
I_avg = mean(abs(err_tx).^2);

SIR_avg_lin = S_avg / I_avg;
SIR_avg_dB = 10*log10(SIR_avg_lin);

fprintf("part a:\n\tavg SIR at TX output = %.4f dB\n", SIR_avg_dB);


% part B
E_s = mean(abs(sym).^2); % average symbol energy
E0 = sum(abs(rrc).^2); % energy in RX filter

% signal power per symbol at MF output (no ISI)
S = abs(gd0)^2 * E_s;
SNR_required_lin = 10^(SNR_required_dB/10);

N_out = S / SNR_required_lin; % desired noise power at MF output
sigma_w2 = N_out / E0; % white noise variance at channel input

% generate complex AWGN at channel
w = sqrt(sigma_w2/2) * (randn(size(x)) + 1j*randn(size(x)));
r = x + w;

% matched filter output
y = filter(rrc, 1, r);

% plot zoomed MF envelope
figure;
plot(abs(y(1:Nplot)), 'LineWidth', 1);
grid on;
title('Q7(b) MF output envelope |y[n]| (first 300 symbols)');
xlabel('Sample index');
ylabel('Envelope');

% average SNIR at MF output
d_tot = span*sps;
idxs_y = d_tot+1 : sps : length(y);
Ns_y = min(length(sym), length(idxs_y));

y_samp = y(idxs_y(1:Ns_y));
sym_use = sym(1:Ns_y);

% ideal RC contribution at MF output
y_ideal = gd0 * sym_use;

% ISI + noise
err_mf = y_samp - y_ideal;

S_avg = mean(abs(y_ideal).^2);
IN_avg = mean(abs(err_mf).^2);
SNIR_lin = S_avg / IN_avg;
SNIR_dB = 10*log10(SNIR_lin);

fprintf("part b:\n\tavg SNIR at MF output = %.4f dB\n", SNIR_dB);


% part C - decode symbols, SER and BER
const = [ ( 1+1j) (-1+1j) (-1-1j) ( 1-1j) ]/sqrt(2);

% inverse gray map corresponding to "const" order above
bits_table = [0 0;   % for const(1)
              0 1;   % for const(2)
              1 1;   % for const(3)
              1 0];  % for const(4)

dec_sym  = zeros(Ns_y,1);
dec_bits = zeros(2*Ns_y,1);

for n = 1:Ns_y
    % nearest neighbor decision
    [~,k] = min(abs(y_samp(n) - const));
    dec_sym(n) = const(k);

    dec_bits(2*n-1:2*n) = bits_table(k,:).';
end

% symbol error rate 
sym_errs = sum(dec_sym ~= sym_use);
SER = sym_errs / Ns_y;

% bit error rate 
bits_use = bits(1:2*Ns_y);
bit_errs = sum(dec_bits ~= bits_use);
BER = bit_errs / (2*Ns_y);

fprintf("part c:\n\t# symbol errors = %d / %d\n", sym_errs, Ns_y);
fprintf("\tSER estimate = %.6e\n", SER);
fprintf("\t# bit errors = %d / %d\n", bit_errs, 2*Ns_y);
fprintf("\tBER estimate = %.6e\n", BER);


% idk why but my MATLAB is not having it with "qfunc" so I just made my own
function [value] = qfunc(x)
    value = 0.5*erfc(x/sqrt(2));
end

% function that computes estimates of P_e as a function of gamma for d_min
% coefficients (a = 2.5, b = 1)
function [p_e] = get_pe_estimate(gamma_k)
    p_e = 2.5 * qfunc(sqrt(gamma_k));
end















