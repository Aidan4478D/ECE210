% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE211 problem_set09.m
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all;

%% Question 1 Part A
N = 1e6;
v = 5; % degrees of freedom
alpha = 0.544;

% gaussian distribution
gaussian = randn(N, 1);

% t-distribution
t = trnd(v, N, 1);
k = 1 / (v / (v - 2)); % scaling factor
t_scaled = t * k; 

% cauchy distribution
U = rand(N, 1);
cauchy = alpha * tan(pi * U);

% fraction of time absolute value < 1
frac_gaussian = sum(abs(gaussian) < 1) / N
frac_t_scaled = sum(abs(t_scaled) < 1) / N
frac_cauchy = sum(abs(cauchy) < 1) / N

figure;
histogram(gaussian);
hold on;
xline(1,'--');
xline(-1, '--');
title('Gaussian Distribution');
xlabel('Value');
ylabel('Probability Density');

figure;
histogram(t_scaled);
hold on;
xline(1, '--');
xline(-1, '--');
title("Scaled t-Distribution");
xlabel('Value');
ylabel('Probability Density');

figure;
histogram(cauchy);
hold on;
xline(1, '--');
xline(-1, '--');
title('Cauchy Distribution');
xlabel('Value');
ylabel('Probability Density');


%% Question 1 Part B
l = 1e5;
seg_n = N / l;
re_gaussian = reshape(gaussian, [l, seg_n]);
re_t_scaled = reshape(t_scaled, [l, seg_n]);
re_cauchy = reshape(cauchy, [l, seg_n]);

u_gaussian = mean(re_gaussian)
u_t_scaled = mean(re_t_scaled)
u_cauchy = mean(re_cauchy)

% due to the heavy tails of the cauchy distribution, the cauchy mean is not
% defined, thus it is not reasonable to say the cauchy mean is equal to 0. 


%% Question 2
