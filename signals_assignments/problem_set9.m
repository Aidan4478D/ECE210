% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE211 problem_set09.m
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all;

%% Question 1 Part A
N = 1e6;
v = 5; % degrees of freedom
alpha = 0.544;

% gaussian distribution
gaussian = randn(N, 1);

% t-distribution
t = trnd(v, N, 1);
k = 1 / sqrt((v / (v - 2))); % scaling factor
t_scaled = t * k; 

% cauchy distribution
U = rand(N, 1);
cauchy = alpha * tan(pi * U);

% fraction of time absolute value < 1
frac_gaussian = sum(abs(gaussian) < 1) / N
frac_t_scaled = sum(abs(t_scaled) < 1) / N
frac_cauchy = sum(abs(cauchy) < 1) / N

figure;
histogram(gaussian);
hold on;
xline(1,'--');
xline(-1, '--');
title('Gaussian Distribution');
xlabel('Value');
ylabel('Probability Density');

figure;
histogram(t_scaled);
hold on;
xline(1, '--');
xline(-1, '--');
title("Scaled t-Distribution");
xlabel('Value');
ylabel('Probability Density');

figure;
histogram(cauchy);
hold on;
xline(1, '--');
xline(-1, '--');
title('Cauchy Distribution');
xlabel('Value');
ylabel('Probability Density');


%% Question 1 Part B
l = 1e5;
seg_n = N / l;
re_gaussian = reshape(gaussian, [l, seg_n]);
re_t_scaled = reshape(t_scaled, [l, seg_n]);
re_cauchy = reshape(cauchy, [l, seg_n]);

u_gaussian = mean(re_gaussian)
u_t_scaled = mean(re_t_scaled)
u_cauchy = mean(re_cauchy)

% due to the heavy tails of the cauchy distribution, the cauchy mean is not
% defined, thus it is not reasonable to say the cauchy mean is equal to 0. 


%% Question 2 Part A

% 1. model is ARMA(2, 2)
% 2. innovations filter
% 3. v[n] = x[n] - 0.4v[n-1] - 0.2v[n-2] - 1.6x[n-1] + 0.81x[n-2]
%
% 4. A(z) = 1 - 1.6z^-1 + 0.81z^-2
%    B(z) = 1 + 0.4z^-1 + 0.2z^-2
%    H(z) = B(z) / A(z)
%         = (1 + 0.4z^-1 + 0.2z^-2) / (1 - 1.6z^-1 + 0.81z^-2)
%
% 5. Sx(w) = var^2 / |H(w)|^2 
%          = 2 * |(1 - 1.6(e^(jw))^-1) + 0.81(e^(jw))^-2) /
%                 (1 + 0.4(e^(jw))^-1) + 0.2(e^(jw))^-2)|^2

% part 6
num = [1, 0.4, 0.2];
den = [1, -1.6, 0.81];

[z,p,k] = tf2zp(num, den);

figure;
zplane(num, den);
title('Zero-poles of H(z)');
% since all zeros and poles are within the unit circle, H is the minimum-phase


%% Question 2 Part B

% a
N = 1e4;
v = sqrt(2) * randn(N, 1);
x = filter(num, den, v);

% b - time averaging to estimate rx
rx = zeros(7,1);
for m = 0:6
    rx(m + 1) = dot(x(1:N-m), x(m+1:N)) / (N - m);
end

% c - stem plot of rx
% since rx(m) = rx(-m)
rx_complete = [flipud(rx(2:end)); rx];

figure;
stem(-6:6, rx_complete);
title('Correlation');

% d - 7x7 Toeplitz matrix R
R_toe = toeplitz(rx)

% e - eigenvalues of R
[eigvec, eigval0] = eig(R_toe);
[eigval, idx] = sort(diag(eigval0),'descend');

pos_def = all(eigval(:)>0)
% since all eigenvalues are > 0, R is positive definite

% f - estimate 6x6 correlation matrix R
X = toeplitz(x, zeros(1, 7));

% correlation matrix R
R_corr = (X' * X) / N


%% Question 2 Part C

% a
[s_est, w] = pwelch(x, hamming(512), 256, 512);

figure;
plot(w,s_est);
title("Estimated PSD");
xlabel('Frequency')
ylabel('PSD')

% b - peak of PSD estimate at frequency w0 in rad
[peak, i] = max(s_est);
w0 = w(i)

% c - pole angle
pole_angle = abs(angle(p))

% the pole angle is super close to w0


%% Question 2 Part D
[a, varv] = aryule(x, 4);

x0 = filter(1, a, v);

rx0 = zeros(7, 1); 
for m = 0:6
    rx0(m + 1) = dot(x0(1:N-m), x0(m+1:N)) / (N - m);
end

rx0_complete = [flipud(rx0(2:end)); rx0];

figure;
stem(-6:6, rx0_complete);
title('Correlation (Part D)');

figure;
hold on;
stem(1:100, x(1:100));
stem(1:100, x0(1:100));
legend('AR Innovations Signal', 'Original signal');
title('Part D Stem Plot of x0 and x');

% plots do not match!

