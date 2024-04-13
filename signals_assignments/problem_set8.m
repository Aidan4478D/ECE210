% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE211 problem_set05.m
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all;


% M x N matrix
% L matrix angle of incidences
% c constant d/lambda
% P power source vector
% Pn noise power vector 
function [s, A] = generate_vectors(M, N, L, c, PdB, PndB)
    
    % steering vector
    s = 1 / sqrt(M) * exp(-1j * 2 * pi * (0:M-1) * c * cos(L));

    % variances for noise and source
    var_s = 10 .^ (PdB / 10);    
    var_n = 10 .^ (PndB / 10);
    
    % generate random signals V and random scaling coefficients B
    B = sqrt(var_s) * randn(M, N) + 1j * randn(M, N) / sqrt(2);
    V = sqrt(var_n) * randn(M, N) + 1j * randn(M, N) / sqrt(2);

    S = zeros(M, N);
    S(:,(0:M)) = s;
   
    A = S .* B + 1/sqrt(M) * V;
end