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

% time index vectors corresponding to data value vector indices
nh = [-1, 0, 1, 2, 3]; 
nx = [-2, -1, 0, 1, 2]; 

h_plot = s_plot(nh, h);
x_plot = s_plot(nx, x);
y_plot = c_plot(h, nh, x, nx);

function c_p = c_plot(x, nx, h, nh)
    % calculate convolution values
    yd = conv(h, x);

    % create discrete time step vector
    yn = (nh(1) + nx(1)):(nh(end) + nx(end));

    c_p = stem(yn, yd); 
    xlim([yn(1) - length(yn) * 0.2, yn(end) + length(yn) * 0.2]);
    ylim([0, max(yd) + max(yd) * 0.25]);
    ylabel('y = x * h');
    xlabel('Discrete time-step n');
    title('Convolution of x(n) and h(n)');
end

% function that creates a stem plot based on 
function s_p = s_plot(n, x)
    s_p = stem(n, x);

    % for aesthetic purposes
    xlim([n(1) - length(n) * 0.2, n(end) + length(n) * 0.2]);
    ylim([0, max(x) + max(x) * 0.25]);
    xlabel('Discrete time-step n');
end



