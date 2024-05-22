% SPDX-License-Identifier: GPL-3.0-or-later
%
% MA224 Probability PS3 -- Question 3.5 
% Copyright (C) 2024 Aidan Cusa <aidan.cusa@cooper.edu>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 

noise_var99 = randn([10000, 1]) * sqrt(0.09);
noise_var999 = randn([10000, 1]) * sqrt(0.061);
noise_var9999 = randn([10000, 1]) * sqrt(0.046);
noise_var99999 = randn([10000, 1]) * sqrt(0.036);
noise_var999999 = randn([10000, 1]) * sqrt(0.031);
noise_var9999999 = randn([10000, 1]) * sqrt(0.026);


figure;
plot(noise_var99)
title("Noise Variance w/ 99% Correction Rate")

figure;
plot(noise_var999)
title("Noise Variance w/ 99.9% Correction Rate")

figure;
plot(noise_var9999)
title("Noise Variance w/ 99.99% Correction Rate")

figure;
plot(noise_var99999)
title("Noise Variance w/ 99.999% Correction Rate")

figure;
plot(noise_var999999)
title("Noise Variance w/ 99.9999% Correction Rate")

figure;
plot(noise_var9999999)
title("Noise Variance w/ 99.99999% Correction Rate")
