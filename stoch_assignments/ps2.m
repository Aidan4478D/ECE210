% Kristof Jablonowski, Aidan Cusa, Jaeho Cho

clc;
clear;

%% Scenario 1
num_samples = 100000;

% continuous uniform distributions for Y and W
Y = unifrnd(-1, 1, [num_samples ,1]);
W = unifrnd(-1, 1, [num_samples ,1]);
X = Y + W; 

% bayes MMSE estimator function (example 8.5)
bayes_estimate = zeros(num_samples, 1);
bayes_estimate(X >= -3 & X < -1) = 1/2 + X(X >= -3 & X < -1)/2;
bayes_estimate(X >= -1 & X < 1) = 0;
bayes_estimate(X >= 1 & X <= 3) = -1/2 + X(X >= 1 & X <= 3)/2;

% linear MMSE estimator (example 8.6)
LMMSE_estimate = (1/5) * X;

% empirical values
empirical_bayes = mean((Y - bayes_estimate).^2);
empirical_LMMSE = mean((Y - LMMSE_estimate).^2);

% theoretical bayes given (pg. 149) and theoretical LMMSE given (pg. 154)
theoretical_bayes = 1/4;
theoretical_LMMSE = 4/15;

% percent error
peBayes = abs(((empirical_bayes - theoretical_bayes) / theoretical_bayes) * 100);
peLMMSE = abs(((empirical_LMMSE - theoretical_LMMSE) / theoretical_LMMSE) * 100);

disp("Bayes MMSE Estimator:");
fprintf("Simulated MSE: %f\n", empirical_bayes);
fprintf("Theoretical MSE: %f\n", theoretical_bayes);
fprintf("Percent Error: %f\n\n", peBayes);

disp("Linear MMSE Estimator:");
fprintf("Simulated MSE: %f\n", empirical_LMMSE);
fprintf("Theoretical MSE: %f\n", theoretical_LMMSE);
fprintf("Percent Error: %f\n\n", peLMMSE);

%% Scenario 2
num_samples = 10000;
mu_Y = 1;

% different variance values for Y and R
var_Y = [0.25, 0.5, 1];
var_R = [0.25, 0.5, 1];

% different numbers of noisy observations
num_noisy_list = [5, 10, 25, 50]; 

% matrices for results
results_empirical = zeros(length(var_Y), length(num_noisy_list));
results_theoretical = zeros(length(var_Y), length(num_noisy_list));

% simulate for each pair of variances and number of observations
for var_idx = 1:length(var_Y)
    for noise_idx = 1:length(num_noisy_list)

        % assuming var_Y and var_R have equal variances
        sigma_Y2 = var_Y(var_idx);
        sigma_R2 = var_R(var_idx);
        num_noisy = num_noisy_list(noise_idx);
        
        % generate samples for Y and noisy samples for X
        Y = normrnd(mu_Y, sqrt(sigma_Y2), [num_samples, 1]);
        X = repmat(Y, 1, num_noisy) + normrnd(0, sqrt(sigma_R2), [num_samples, num_noisy]);
        
        % linear estimator (Y hat) and empirical MSE
        estimate = (sigma_R2 * mu_Y + sigma_Y2 * sum(X, 2)) / (num_noisy * sigma_Y2 + sigma_R2);
        empirical_MSE = mean((Y - estimate).^2);
        
        % empirical and theoretical MSE (pg. 160)
        theoretical_MSE = sigma_Y2 * sigma_R2 / (num_noisy * sigma_Y2 + sigma_R2);
        
        results_empirical(var_idx, noise_idx) = empirical_MSE;
        results_theoretical(var_idx, noise_idx) = theoretical_MSE;
    end
end

figure; 
hold on;

colors = lines(length(var_Y));
for i = 1:length(var_Y)
    plot(num_noisy_list, results_empirical(i, :), 'o-', 'Color', colors(i,:), 'DisplayName', "Empirical, var_Y = " + var_Y(i) + ", var_R = " + var_R(i));
    plot(num_noisy_list, results_theoretical(i, :), 'x--', 'Color', colors(i,:), 'DisplayName', "Theoretical, var_Y = " + var_Y(i) + ", var_R = " + var_R(i));
end

xlabel('Num Observations');
ylabel('MSE');
title('Empirical vs Theoretical MSE for Noisy Observations');
legend('show');
grid on;


%% Scenario 3
data = load('SATs.mat');
SAT_Verbal = data.SAT_Verbal;
SAT_Math = data.SAT_Math;

% Remove NaN values from the data
validIndices = ~isnan(SAT_Math) & ~isnan(SAT_Verbal);
SAT_Math = SAT_Math(validIndices);
SAT_Verbal = SAT_Verbal(validIndices);

% total SAT score without empty rows
SAT_Total = SAT_Math + SAT_Verbal;
all_scores = true(size(SAT_Total));
mid_scores = (SAT_Total >= 1150 & SAT_Total <= 1250);
high_scores = SAT_Total > 1320;

% parameters for each subset
[slope_all, intercept_all] = calculate_line(SAT_Math(all_scores), SAT_Verbal(all_scores));
[slope_mid, intercept_mid] = calculate_line(SAT_Math(mid_scores), SAT_Verbal(mid_scores));
[slope_high, intercept_high] = calculate_line(SAT_Math(high_scores), SAT_Verbal(high_scores));

% line for the entire range of SAT Math scores
x_line = linspace(min(SAT_Math), max(SAT_Math), 1000);

y_line_all = slope_all * x_line + intercept_all;
y_line_mid = slope_mid * x_line + intercept_mid;
y_line_high = slope_high * x_line + intercept_high;

% plot scores and lines
figure;
hold on;

scatter(SAT_Math(all_scores), SAT_Verbal(all_scores), 'filled', 'DisplayName', 'All Data');
scatter(SAT_Math(mid_scores), SAT_Verbal(mid_scores), 'filled', 'DisplayName', 'Scores 1150-1250');
scatter(SAT_Math(high_scores), SAT_Verbal(high_scores), 'filled', 'DisplayName', 'Scores >1320');

plot(x_line, y_line_all, 'k-', 'DisplayName', 'Least Fit Line - All Scores');
plot(x_line, y_line_mid, 'g-', 'DisplayName', 'Least Fit Line - 1150-1250');
plot(x_line, y_line_high, 'r-', 'DisplayName', 'Least Fit Line - >1320');

xlabel('SAT Math Scores');
ylabel('SAT Verbal Scores');
title('Least Fit Lines of SAT Verbal Scores');
legend('show');
grid on;
hold off;

% when plotting over all data, the estimator reflects a positive
% relationship between math and verbal scores (students who score higher in
% math tend to score higher in verbal) 

% on specific ranges the estimator reflects a negative relationship between 
% math and verbal where an increase in math implies a decrease in verbal
% score. Thus when the overall score is held nearly constant, a higher
% performance in one section may come at the cost of the other

% calculate slope, intercept, and correlation
function [slope, intercept] = calculate_line(x, y)
    % calculate pearson correlation coefficient
    p = corrcoef(x, y);
    % get pearson correlation coefficient between x (math) and y (verbal)
    corr = p(1, 2);

    % slope = cov(x, y) / var(x)
    % since p = cov(x, y) / (std(x)std(y))
    % slope = p * std(y) / (std(x)
    slope = corr * std(y) / std(x);
    intercept = mean(y) - slope * mean(x);
end































