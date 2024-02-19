% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE211 Signal Processing - Problem Set 3
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 

% data values vectors
h = [2, 3, 4, 2, 5];
x = [3, 4, 5, 1, 2]; 

% time index vector // first entry represents n = 0
nh = [-1, 0, 1, 2, 3]; 
nx = [-2, -1, 0, 1, 2]; 

% calculate convolution values
yd = conv(h, x);
% create discrete time step vector
yn = (nh(1) + nx(1)):(nh(end) + nx(end));

% create a figure and turn hold off
figure;
hold off; 
xlabel('Discrete time-step n');

% create stem plot for h(n)
subplot(3, 1, 1);
stem(nh, h);
ylabel('h');
xlim([nh(1) - 2, nh(end) + 2]);
ylim([0, max(h) + max(h) * 0.25]);
title('Plot of h(n)')

% create stem plot for h(n)
subplot(3, 1, 2);
stem(nx, x);
ylabel('x');
xlim([nx(1) - 2, nx(end) + 2]);
ylim([0, max(x) + max(x) * 0.25]);
title('Plot of x(n)')

% create stem plot for convolution y(n) = h * x
subplot(3, 1, 3);
stem(yn, yd)
xlim([yn(1) - 2, yn(end) + 2]);
ylim([0, max(yd) + max(yd) * 0.25]);
ylabel('y = h * x');
title('Convolution of x(n) and h(n)');




