% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE300 Communication Theory - Problem Set 5
% Copyright (C) 2025 Aidan Cusa <aidancusa@gmail.com>

clc; clear; close all; 

% Q1

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


% Q2

% setting E0 = 1 for ease of calculations
constellation = [-3 - 1j, -3 + 1j, ...
                 -1 - 1j, -1 + 1j, ...
                  3 - 1j, 3 + 1j, ...
                  1 - 1j, 1 + 1j];

% part B

display(constellation)

E_s = mean(abs(constellation).^2);
E_b = E_s / 3; % 3 bits / sym

fprintf("E_b = %f * E_0\n", E_b);
%display(E_b)

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
fprintf('total pairs: %d (should 56)\n\n', sum(abs(distances) > 0));


% idk why but my MATLAB is not having it with "qfunc" so I just made my own
function [value] = qfunc(x)
    value = 0.5*erfc(x/sqrt(2));
end















