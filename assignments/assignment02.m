% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment02.m -- A New Way of Thinking
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current workspace

%% part 1
u = -4:2:4;
v = 0:(pi/4):pi;

%% part 2
f = prod(1:10);

%% part 3

% part A
A = zeros(2, 4);
A(1, 1) = 1;
A(2, 3) = 1; 

% display(A)

% part B
b = reshape(1:16, 2, [])';
B = reshape(b, [4, 4]);

%% part 4

t = linspace(-pi, pi, 50);



