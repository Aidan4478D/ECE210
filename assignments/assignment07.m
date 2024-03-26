% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment08.m -- Under Pressure
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 

fs = 80e3; 

%% part 1
brake_data = readmatrix("40p_1000ms.csv");
%norm_brake_data = rescale(brake_data(:, 2))%normc(brake_data)

%% part 2

t = normalize(brake_data(:, 1), 'range', [0, 1]);
norm_brake_data = normalize(brake_data(:, 2), 'range', [0, 1]);

x = fft(norm_brake_data);
f = fs * linspace(0, pi, numel(x)) / (2 * pi); % normalize fs to w range


db = 20 * log10(abs(x)); % get db value

figure;
plot(f, db); 