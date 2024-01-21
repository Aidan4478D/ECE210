% SPDX-License-Identifier: GPL-3.0-or-later
%
% Copyright (C) 2024 Aidan Cusa <aaidancusa@gmail.com>


clc;
clear;

%% part 1
u = [11, 13, 17];
v = [-1; -1; -1];
A = [-1 * u; 2 * u; 7 * u];
B = [A', v];

%% part 2
c = exp(1j * pi);
d = sqrt(1j);
l = floor((8.4018 * 10^6)^(2.1 / 2));
k = floor(100 * log(2)) + ceil(exp(7.5858));

%% part 3
A1 = [1, -11, -3
      1,   1,  0
      2,   5   1];

b = [-37; -1; 10];

% method #1 of mldivide
% x = A1\b 

% method #2 of mldivide
x = mldivide(A1, b);
