% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE211 problem_set05.m
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all;

%% Part 1
M = 100;
N = 200;
L = 3;
c = 0.5; % d/lambda

PdB = [0, -2, -4];
PndB = 10;
thetas1 = [10, 25, 70] * 2 * pi / 180;
thetas2 = [10, 12, 70] * 2 * pi / 180;

% data matrices A
A1 = generate_data(M, N, thetas1, c, PdB, PndB);
A2 = generate_data(M, N, thetas2, c, PdB, PndB);

% correlation matrices 
R1 = A1 * A1' / N;
R2 = A2 * A2' / N;

%% part 2
[svals1, eigvals1, U1] = SVD_EDecomp(A1, R1);
[svals2, eigvals2, U2] = SVD_EDecomp(A2, R2);

% plot singular values
figure;
subplot(2, 1, 1);
stem(svals1);
title("Singular Values Stem plot A1");

subplot(2, 1, 2);
stem(svals2);
title("Singular Values Stem plot A2");

% plot sorted eigenvalues of R
figure;
subplot(2, 1, 1);
stem(eigvals1);
title("Eigenvalues Stem plot R1");

subplot(2, 1, 2);
stem(eigvals2);
title("Eigenvalues Stem plot R2");

sval_ratio1 = svals1(3)/svals1(4)
sval_ratio2 = svals2(3)/svals2(4) 

eigval_ratio1 = eigvals1(3)/eigvals1(4)
eigval_ratio2 = eigvals2(3)/eigvals2(4)




% M x N matrix
% L matrix angle of incidences
% c constant d/lambda
% P power source vector
% Pn noise power vector 
function A = generate_data(M, N, theta, c, PdB, PndB)
    
    % num sources
    L = length(theta);
    
    % steering vector
    s = exp(-1j * 2 * pi * (0:M-1)' * (c * cos(theta)));
    s = s / sqrt(M);  % Normalize
    
    % variances
    var_s = 10 .^ (PdB / 10);
    var_n = 10 ^ (PndB / 10);
    
    B = sqrt(var_s') .* randn(L, N) + 1j .* randn(L, N) / sqrt(2);
    V = sqrt(var_n) .* randn(M, N) + 1j .* randn(M, N) / sqrt(2);
    
    % Add noise, scaled by sqrt of M
    A = s * B + V / sqrt(M);
end


function [sval, eigval, U] = SVD_EDecomp(A, R)
    [U, sval, V] = svd(A);
    [eigvec, eigval_] = eig(R);
    [eigval, idx] = sort(diag(eigval_), 'descend');
    eigvec = eigvec(:,idx);
end
