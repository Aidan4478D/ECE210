% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment01.m -- Practice and Perform Basic MATLAB Operations
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current workspace

%% part 1
u = [11, 13, 17];            % create row vector
v = [-1; -1; -1];            % create column vector

A0 = [-1 * u; 2 * u; 7 * u]; % create matrix by multiplying existing row
                             % vector u by a scalar

B = [A0', v];                % create a matrix using the transpose of
                             % matrix A and column vector v

%% part 2
c = exp(1j * pi / 4);
d = sqrt(1j);
l = nthroot(8.4018 * 10^6, 2.1);
k = floor(100 * log(2)) + ceil(exp(7.5858));

%% part 3
A = [1, -11, -3
     1,   1,  0
     2,   5   1];

b = [-37; -1; 10];

% going to use mldivide which uses matrix left division in order to solve
% systems of linear equations of the form Ax = b

% method #1 of mldivide
% x = A\b

% method #2 of mldivide
x = mldivide(A, b);


OUTPUT = 'output.pdf';

r = groot;

delete(OUTPUT);
for i = numel(r.Children):-1:1
    exportgraphics(r.Children(i), OUTPUT, 'Append', true, 'ContentType', 'vector');
end
