% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment02.m -- A New Way of Thinking
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 

%% part 1
u = -4:2:4;
v = 0:(pi/4):pi; % could also use linspace(0, pi, 5) but meh

%% part 2
f = prod(1:10);

%% part 3

% part A
A = zeros(2, 4);
A(1, 1) = 1;
A(2, 3) = 1; 

% part B
b = reshape(1:16, 2, [])'; % reshape 1-16 into two rows then transpose it to 
                           % become one even column and one odd column

B = reshape(b, [4, 4]);    % now reshape the two columns into a 4x4 matrix
                           % will get desired result due to column major
                           % ordering

%% part 4
n = 0:50; 
t = linspace(-pi, pi, 1000);

a_n = (2 * n + 1)';   % transpose in order to make a_n into a column vector
                      % so that it can be multiplied by t properly

s = sum(sin(a_n * t) ./ a_n);  % element wise division for dimensions to
                               % match up

plot(t, s)

