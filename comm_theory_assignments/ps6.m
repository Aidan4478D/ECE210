% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE300 Communication Theory - Problem Set 5
% Copyright (C) 2025 Aidan Cusa <aidancusa@gmail.com>

clc; clear; close all; 

%% Q1
fprintf("question 1\n")

% need to find k to satisfy Hamming bound
n = 15;
t = 3;

% sum of binomial coefficients
binom_sum = 0;
for i = 0:t
    binom_sum = binom_sum + nchoosek(n, i);
end
fprintf("\tneed to find k such that 2^(n-k) >= %d\n", binom_sum);

% find max k that satisfies Hamming bound
max_k = 0;
for k = 1:n
    val = 2^(n-k);
    if binom_sum <= val
        max_k = k;
    end
end
fprintf('\tlargest k satisfying Hamming bound: %d\n', max_k);

%% Q2
fprintf("question 2\n")

% equation on deck 3, slide 5
n = 20;
l = 6;

num_burst_errors = (n - l + 1) * 2^(l - 2);
fprintf('\tnumber of burst errors: %d\n', num_burst_errors);

%% Q3
fprintf("question 3\n")

G = [
    1 0 0 1 1 0 1 0 1;
    0 1 0 0 1 1 1 0 1;
    0 0 1 0 1 0 0 1 0
];

% k = rows, n = columns
[k, n] = size(G);

% part A:
fprintf("\t(a)k=%d\n", k)
fprintf("\tn=%d\n", n)

% part B:
% should be 2^3=8
num_codewords = 2^k;
fprintf("\t(b) num codewords: %d\n", num_codewords)

% part C:
% codewords calculated by hand
codewords = [
    0 0 0 0 0 0 0 0 0;
    0 0 1 0 1 0 0 1 0;
    0 1 0 0 1 1 1 0 1;
    0 1 1 0 0 1 1 1 1;
    1 0 0 1 1 0 1 0 1;
    1 0 1 1 0 0 1 1 1;
    1 1 0 1 0 1 0 0 0;
    1 1 1 1 1 1 0 1 0;
];
%disp(codewords)

% find d_min
weights = sum(codewords, 2);
d_min = min(weights(weights > 0)); % ignore the all-zero codeword (first row)
t = floor((d_min - 1) / 2); % deck 1 side 19

fprintf('\t(c) d_min = %d\n', d_min);
fprintf('\terror correcting capability t = %d\n', t);

% part D:
% from deck 2 slide 46
P = G(:, 4:end);
H = [P', eye(n-k)]; 
fprintf('\t(d) parity check matrix H:');

% should be all zeros
disp(H);
fprintf('\tcheck GH^T = 0: %d\n', sum(sum(mod(G * H.', 2))))  

% part E:
d_vec = [1 1 0];
c_vec = mod(d_vec * G, 2);
fprintf('\t(e) encoded vector for [1 1 0]:');
disp(c_vec);

% part F:
fprintf('\t(f)original codeword c:'); 
disp(c_vec);

% errors in first and last position
e_true = [1 0 0 0 0 0 0 0 1];
r = mod(c_vec + e_true, 2);
fprintf('\treceived vector r:'); 
disp(r);

% syndrome
s = mod(r * H.', 2);
fprintf('\tsyndrome s:');
disp(s);

% enumerate all codewords
codewords = zeros(num_codewords, 9);
d_list = zeros(num_codewords, 3);

idx = 1;
for d1 = 0:1
    for d2 = 0:1
        for d3 = 0:1
            dd = [d1 d2 d3];
            d_list(idx,:) = dd;
            codewords(idx,:) = mod(dd * G, 2);
            idx = idx + 1;
        end
    end
end

fprintf('\tall codewords (u | c):\n');
disp([d_list codewords]);

% each row is r + c_i
coset = mod(r + codewords, 2);

% Hamming weight of each coset vector
weights = sum(coset, 2);

fprintf('\n\tcoset entries (x = r + c_i) and their weights:\n');
for i = 1:num_codewords
    fprintf('\t\tc_%d = ', i-1);
    disp(codewords(i,:));
    fprintf('\t\tx_%d = ', i-1);
    disp(coset(i,:));
    fprintf('\t\tweight(x_%d) = %d\n', i-1, weights(i));
end

% identify coset leader (minimum-weight vector in coset)
[min_w, min_idx] = min(weights);
e_leader = coset(min_idx, :);

fprintf('\n\tcoset leader:\n');
fprintf('\tindex in coset: %d\n', min_idx);
fprintf('\tcoset leader e_leader:'); 
disp(e_leader);
fprintf('\tweight(e_leader) = %d\n', min_w);

% decode
c_hat = mod(r + e_leader, 2);
fprintf('\n\tdecoded codeword:'); disp(c_hat);

% check if we recovered the original codeword
is_correct = isequal(c_hat, c_vec);
fprintf('\tdid we recover the original codeword: %d\n', is_correct);


%% Q4

% part A: 
% code rate = bits clocked in / bits cloced out = 2/3
% constrain length (L) = max time depth in system = 5
% number of states = 2^((L-1) * k) = 256

% part B:
% deck 5 slide 19
% "for a rate k/n code, for each present state, there are exactly 2^k
% previous states"
% thus, there are 2^2 = 4 previous states

% part C:
top = [0 1 1 0]; % [T1, T2, T3, T4]
bottom = [1 1 1 0]; % [B1, B2, B3, B4]

% bits T2 T3 and T4 of the present top state were the bits T1 T2 and T3 of
% the previous state ==> previous top must be [1 1 0 X] and the previous 
% bottom must be [1 1 1 Y]

% the inputs that caused this transition are the values currently in the 
% first position -> Input top = 0, input bottom = 1

% thus, the previous possible states are:
% top: [1 1 0 0] bottom: [1 1 1 0]
% top: [1 1 0 0] bottom: [1 1 1 1]
% top: [1 1 0 1] bottom: [1 1 1 0]
% top: [1 1 0 1] bottom: [1 1 1 1]


















