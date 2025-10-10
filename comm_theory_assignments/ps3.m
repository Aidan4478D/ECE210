% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE300 Communication Theory - Problem Set 2
% Copyright (C) 2025 Aidan Cusa <aidancusa@gmail.com>

clc; clear; close all; 

%% Q1
disp("Question 1:")

% part a
probs = [0.225, 0.205, 0.19, 0.14, 0.08, 0.07, 0.05, 0.04];

code_lengths = ceil(-log2(probs));
disp("Code Lengths")
disp(code_lengths)

kraft_sum = sum(2.^(-code_lengths));
fprintf("Kraft sum: %.3f\n", kraft_sum)

kraft_satisfied = (kraft_sum <= 1);
fprintf("Kraft inequality satisfied: %d\n", kraft_satisfied)

% part b
entropy = sum(-probs .* log2(probs));
avg_lengths = sum(probs .* code_lengths);

fprintf("Entropy of source H(X): %.3f\n", entropy);
fprintf("Average code lengths: %.3f\n", avg_lengths)

inequality = ((entropy <= avg_lengths) < entropy + 1);
fprintf("Verification of inequality: %d\n", inequality)

%% Q2
disp("Question 2:")

% code lengths in tree are smaller than in the previous problem!
huffman_lengths = [2, 2, 3, 3, 4, 4, 4, 4];
kraft_sum = sum(2.^(-huffman_lengths));
fprintf("Kraft sum (huffman): %.3f\n", kraft_sum)

kraft_satisfied = (kraft_sum <= 1);
fprintf("Kraft inequality satisfied (huffman): %d\n", kraft_satisfied)

% part b
entropy = sum(-probs .* log2(probs));
avg_lengths = sum(probs .* huffman_lengths);

fprintf("Entropy of source H(X) (huffman): %.3f\n", entropy);
fprintf("Average code lengths (huffman): %.3f\n", avg_lengths)

inequality = ((entropy <= avg_lengths) < entropy + 1);
fprintf("Verification of inequality (huffman): %d\n", inequality)

%% Q3 - this took so long
disp("Question 3:")

% create dictionary for easier calculation
codes = ["000", "001", "010", "011", "100", "101", "110", "111"];
probs = [0.08815, 0.11685, 0.09675, 0.12825, 0.09, 0.14, 0.15, 0.19]; % calculated by hand
d = dictionary(codes, probs);

% sum of probabilities should be one
% check = sum(probs); 
% disp(check)

% P(Y3)
p_y3_0 = d("000") + d("010") + d("100") + d("110");
p_y3_1 = d("001") + d("011") + d("101") + d("111");

q3_1_results = [p_y3_0; p_y3_1];
q3_1_labels = ['P(Y3 = 0)'; 'P(Y3 = 1)'];
results_table_1 = table(q3_1_labels, q3_1_results, 'VariableNames', {'Event','Probability'});
disp(results_table_1);

entropy_y3 = -(p_y3_0 * log2(p_y3_0) + p_y3_1 * log2(p_y3_1));
fprintf("Entropy H(Y3): %.5f\n", entropy_y3)


% P(Y3 | Y2)
% marginals
p_y2_0 = d("000") + d("001") + d("100") + d("101");
p_y2_1 = d("010") + d("011") + d("110") + d("111");

% conditionals
p_y3y2_00 = (d("000") + d("100")) / p_y2_0; % P(Y3=0 | Y2=0)
p_y3y2_10 = (d("001") + d("101")) / p_y2_0; % P(Y3=1 | Y2=0)
p_y3y2_01 = (d("010") + d("110")) / p_y2_1; % P(Y3=0 | Y2=1)
p_y3y2_11 = (d("011") + d("111")) / p_y2_1; % P(Y3=1 | Y2=1)

entropy_y3y2_0 = -(p_y3y2_00 * log2(p_y3y2_00) + p_y3y2_10 * log2(p_y3y2_10));
entropy_y3y2_1 = -(p_y3y2_01 * log2(p_y3y2_01) + p_y3y2_11 * log2(p_y3y2_11));

q3_2_results = [p_y3y2_00; p_y3y2_01; p_y3y2_10; p_y3y2_11];
q3_2_labels = ['P(Y3=0|Y2=0)'; 'P(Y3=0|Y2=1)'; 'P(Y3=1|Y2=0)'; 'P(Y3=1|Y2=1)'];
results_table_1 = table(q3_2_labels, q3_2_results, 'VariableNames', {'Event','Probability'});
disp(results_table_1);

entropy_y3y2 = (p_y2_0 * entropy_y3y2_0 + p_y2_1*entropy_y3y2_1);
fprintf("Entropy H(Y3|Y2): %.5f\n", entropy_y3y2)

% P(Y3 | Y1 Y2)
% marginals
p_y1y2_00 = d("000") + d("001");
p_y1y2_01 = d("010") + d("011");
p_y1y2_10 = d("100") + d("101");
p_y1y2_11 = d("110") + d("111");

% P(Y3 | Y1=0, Y2=0)
p_y3y1y2_000 = d("000") / p_y1y2_00;
p_y3y1y2_100 = d("001") / p_y1y2_00;
entropy_y3y1y2_00 = -(p_y3y1y2_000*log2(p_y3y1y2_000) + p_y3y1y2_100*log2(p_y3y1y2_100));

% P(Y3 | Y1=0, Y2=1)
p_y3y1y2_001 = d("010") / p_y1y2_01;
p_y3y1y2_101 = d("011") / p_y1y2_01;
entropy_y3y1y2_01 = -(p_y3y1y2_001*log2(p_y3y1y2_001) + p_y3y1y2_101*log2(p_y3y1y2_101));

% P(Y3 | Y1=1, Y2=0)
p_y3y1y2_010 = d("100") / p_y1y2_10;
p_y3y1y2_110 = d("101") / p_y1y2_10;
entropy_y3y1y2_10 = -(p_y3y1y2_010*log2(p_y3y1y2_010) + p_y3y1y2_110*log2(p_y3y1y2_110));

% P(Y3 | Y1=1, Y2=1)
p_y3y1y2_011 = d("110") / p_y1y2_11;
p_y3y1y2_111 = d("111") / p_y1y2_11;
entropy_y3y1y2_11 = -(p_y3y1y2_011*log2(p_y3y1y2_011) + p_y3y1y2_111*log2(p_y3y1y2_111));

q3_3_results = [p_y3y1y2_000; p_y3y1y2_001; p_y3y1y2_010; p_y3y1y2_011; p_y3y1y2_100; p_y3y1y2_101; p_y3y1y2_110; p_y3y1y2_111];
q3_3_labels = ['P(Y3=0|Y1=0,Y2=0)'; 'P(Y3=0|Y1=0,Y2=1)'; 'P(Y3=0|Y1=1,Y2=0)'; 'P(Y3=0|Y1=1,Y2=1)'; ...
               'P(Y3=1|Y1=0,Y2=0)'; 'P(Y3=1|Y1=0,Y2=1)'; 'P(Y3=1|Y1=1,Y2=0)'; 'P(Y3=1|Y1=1,Y2=1)'];

results_table_1 = table(q3_3_labels, q3_3_results, 'VariableNames', {'Event','Probability'});
disp(results_table_1);

entropy_y3y1y2 = p_y1y2_00 * entropy_y3y1y2_00 + ...
                 p_y1y2_01 * entropy_y3y1y2_01 + ...
                 p_y1y2_10 * entropy_y3y1y2_10 + ...
                 p_y1y2_11 * entropy_y3y1y2_11;

fprintf("Entropy H(Y3|Y2): %.5f\n", entropy_y3y1y2)


%% Q4
disp("Question 4:")

sig_squared = [6, 2, 1.6, 0.9];
target_R = 2;

lamb_low  = 0;
lamb_high = max(sig_squared);  % upper bound for lambda

compute_R = @(lamb) sum( 0.5 * log2( sig_squared(sig_squared > lamb) ./ lamb ) );

stop_lambda_tol = 1e-5;
stop_rate_tol = 1e-7;

while ( (lamb_high - lamb_low) > stop_lambda_tol )
    lamb = 0.5*(lamb_low + lamb_high);
    R = compute_R(lamb);
    
    if abs(R - target_R) <= stop_rate_tol % good enough
        lamb_low = lamb; 
        lamb_high = lamb;
    elseif R > target_R % rate too big -> increase lambda
        lamb_low = lamb;
    else
        lamb_high = lamb;  % rate too small -> decrease lambda
    end
end

lamb_final = 0.5*(lamb_low + lamb_high);

% distortions and rates per component
D_i = min(sig_squared, lamb_final);
R_i = 0.5 * log2(sig_squared ./ lamb_final) .* (sig_squared > lamb_final);

D = sum(D_i);
R = sum(R_i);

fprintf('Optimal lambda: %.6f\n', lamb_final);
fprintf('Rate R = %.6f bits,  Distortion D = %.6f\n', R, D);

% per-component allocations
fprintf('\nDistortion and rate per component:\n');
fprintf(' i     sigma^2    D_i        R_i\n');
for i = 1:numel(sig_squared)
    fprintf('%2d    %7.4f    %7.4f    %7.4f\n', i, sig_squared(i), D_i(i), R_i(i));
end

%% Q5
disp("Question 5:")

% Using Shannon Formula:

N_0 = 1e-8; % 10e-8 Joules
B = 1e6; % 1 MHz
C = 1e6; % 1 Mbps

P = (2^(C / B) - 1) * N_0 * B;
fprintf("Transmission power (W): %0.4f\n", P)
