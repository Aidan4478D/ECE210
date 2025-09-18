% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE300 Communication Theory - Problem Set 1
% Copyright (C) 2025 Aidan Cusa <aidancusa@gmail.com>

clc; clear; close all; 

% 8-PSK
k = 0:7;
theta = 2*pi*k/8;
A = 1 / sqrt(2 - 2 * cos(pi / 4)); % calculated in part 3
S_8psk = A * [cos(theta); sin(theta)];   % 2 x 8

% 16-QAM
lev = [-1.5 -0.5 0.5 1.5];
[X, Y] = meshgrid(lev, lev);
S_16qam = [X(:)'; Y(:)'];   % 2 x 16

% 4-Orthogonal (N=4, M=4)
A = 1/sqrt(2);
S_4orth = A * eye(4);

% Run your function and compute dB
constellations = {'8-PSK', '16-QAM', '4-Orthogonal'};
[eb_8psk, se_8psk] = const_metrics(S_8psk);
[eb_16qam, se_16qam] = const_metrics(S_16qam);
[eb_4o, se_4o] = const_metrics(S_4orth);

eb_list = [eb_8psk; eb_16qam; eb_4o];
se_list = [se_8psk; se_16qam; se_4o];

% results in a table
results_table = table(constellations', eb_list, se_list, 'VariableNames', {'Constellation', 'Eb_dB_per_bit', 'Spectral_Efficiency_bits_dim'});
disp(results_table);

% part d
%   16-QAM is the most spectrally efficient (2 bits/dim)
%   8-PSK is the next most spectrally efficient (1.5 bits/dim)
%   4-orthogonal is the least spectrally efficient (0.5 bits/dim)

% part e
%   most power efficient is 16-QAM (-2.0412)
%   difference btwee nmost and least is 16-QAM - 4-orthogonal 
%       ==> (-2.0412 - (-6.0206)) = 3.9794 dB/bit


% S is NxM, cols are constellation points in R^N
function [energy_bit_dB, spectral_efficiency] = const_metrics(S)

    [N, M] = size(S);

    % average energy per symbol
    energy_symbol = mean(sum(S.^2, 1));

    % energy per bit
    energy_bit = 1 / log2(M) * energy_symbol;
    energy_bit_dB = 10 * log10(energy_bit);
    
    % spectral efficiency
    spectral_efficiency = log2(M) / N;
end