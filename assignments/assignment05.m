% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment05.m -- Plotting, Scheming Even
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 

%% part 1
n = 50; 
t = linspace(-pi, pi, 1000);
                             
figure; 
hold on;

for i = 0:n
    plot(t, fourier_sum(0:i, t))
end

title('Square Wave with N=50');
xlabel('Time [t]');
xlim([-pi, pi])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
xticks([-pi, -pi/2, 0, pi/2, pi]);
ylabel('y');


function f_s = fourier_sum(n, t)
    a_n = (2 * n + 1)';
    f_s = sum(sin(a_n * t) ./ a_n);
end

%% part 2