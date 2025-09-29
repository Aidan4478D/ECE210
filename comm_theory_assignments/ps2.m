% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE300 Communication Theory - Problem Set 2
% Copyright (C) 2025 Aidan Cusa <aidancusa@gmail.com>

clc; clear; close all; 

%% Part A

% A determinstic decision rule is just a function from two possible
% observations (Red2, Blue2) to two possible actions (Red1, Blue1) so it's
% fully represented by a 2x1 matrix.

% The representation of this matrix is within the "MAP_decision" and 
% "ML_decision" functions.

%% Part D 
b1_cases = [2, 6, 2, 6];
r1_cases = [8, 4, 8, 4];
b2_cases = [2, 2, 6, 6];
r2_cases = [8, 8, 4, 4];

% DECISION RULES:
% [1,1]: Always choose red 1
% [1,2]: Match what you see (red2 -> red1 and blue2 -> blue1)
% [2,1]: Flip what you see (red2 -> blue1 and blue2 -> red1)
% [2,2]: Always choose blue 1


% Comparison between MAP and ML
% Case 1: (r1=8,b1=2,r2=8,b2=2)
% MAP has decision rule [1, 1], error 0.200 (expimental 0.2013)
% ML has decision rule [1, 2], error 0.2909 (expimental 0.2917)
% Here, the prior favors red and urn 2 is also read leaning but ML flips to
% blue1 when blue2 is seen causing extra erros because red1 is so frequent.
% In this case, MAP "trusts" the prior on both branches and results in less
% erorr. 

% Case 2: (r1=4,b1=6,r2=8,b2=2)
% MAP has decision rule [2, 2], error 0.400 (expimental 0.4017)
% ML has decision rule [1, 2], error 0.5091 (expimental 0.5093)
% Here, the prior favors blue but urn 2 is red leaning. Because of this,
% ML chooses red1 on red2 a lot which leads to more error as blue1 is
% actually more common. MAP sticks to the prior and avoids those flips
% which leads to less error.

% Case 3: (r1=8,b1=2, r2=4,b2=6)
% MAP has decision rule [1, 1], error 0.200 (expimental 0.2003)
% ML has decision rule [1, 2], error 0.5091 (expimental 0.5066)
% Here, there is a large gap as the prior leans heavily on red but urn 2
% now leans mroe blue. ML picks blue1 on blue2 and gets it wrong most often
% as most first draws were red1. MAP will always pick red1 which makes its
% error equal to the prior which destroys ML here (lol). 

% Case 4: (r1=4,b1=6,r2=4,b2=6)
% MAP has decision rule [2, 2], error 0.400 (expimental 0.4014)
% ML has decision rule [1, 2], error 0.4364 (expimental 0.4365)
% Herem the prior favors blue and urn 2 also leans towards blue. ML still
% flips to red1 on red2 while MAP stays blue1 on both. The gap is smaller
% than case 2 because the likelihood isn't fighting the prior as hard, but
% MAP is still better.

% In general with these four setups, MAP has a lower error rate than ML
% every time which makes it a better decision rule. MAP is significantly 
% better than ML when there is a strong prior probability favoring one
% class which leads MAP to almost always choose the more probable class. 
% MAP is only slightly better than ML when the likelihoods are identical 
% or the difference in error comes from a single decision (ex. case 1 where
% the decision for blue results in a different error.), indicating the 
% priors are only slightly unequal or the likelihoods for the two classes 
% are similar. However, if the priors were equal (or if both rules choose 
% the same action on each observation) then MAP would be the same as ML. 
% This equality is not borne out in any of the four cases, confirming that 
% the MAP calculations used unequal priors. Finally, the experimental 
% resutls all conformed with the "theoretical" probabilties of error as 
% seen above. 


%% Part B
% just going to assume that r1, b1, r2, and b2 all have the same length
for case_idx = 1:length(r1_cases)
    r1 = r1_cases(case_idx);
    b1 = b1_cases(case_idx);
    r2 = r2_cases(case_idx);
    b2 = b2_cases(case_idx);
    
    fprintf('\n\nCase %d: b1=%d, r1=%d, b2=%d, r2=%d\n', case_idx, b1, r1, b2, r2);
    
    % likelihoods of second draw given the first

    % +1 as if urn 1 was red1 you would move a red ball into urn 2 before
    % drawing so urn 2 now has r2 + 1 reds and similar reasoning for blue
    likelihoods = [(r2 + 1), b2; r2, (b2 + 1)] / (r2 + b2 + 1);

    % priors over urn 1 color
    P_prior = [r1, b1, 0, 0] / (r1 + b1);

    % marginals of the second draw P(Red2) and P(Blue2)
    P_prior(3) = likelihoods(1,1) * P_prior(1) + likelihoods(2,1) * P_prior(2);
    P_prior(4) = likelihoods(1,2) * P_prior(1) + likelihoods(2,2) * P_prior(2);
    
    % posteriors over urn 1 given what you saw in urn 2 (bayes theorem)
    Posterior = [likelihoods(1,1) * P_prior(1) / P_prior(3), likelihoods(2,1) * P_prior(2) / P_prior(3); ...
                 likelihoods(1,2) * P_prior(1) / P_prior(4), likelihoods(2,2) * P_prior(2) / P_prior(4)];
    
    % calculate MAP and ML decision rules
    [map_decision_rule] = MAP_decision(Posterior);
    [ml_decision_rule] = ML_decision(likelihoods);
    
    % calculate errors
    [map_error_Red2, map_error_Blue2] = error_prob(map_decision_rule, P_prior, likelihoods);
    [ml_error_Red2, ml_error_Blue2] = error_prob(ml_decision_rule, P_prior, likelihoods);
    
    
    fprintf('\nMAP Error Probability (Red2 observed): %f\n', map_error_Red2);
    fprintf('MAP Error Probability (Blue2 observed): %f\n', map_error_Blue2);
    fprintf('\nML Error Probability (Red2 observed): %f\n', ml_error_Red2);
    fprintf('ML Error Probability (Blue2 observed): %f\n', ml_error_Blue2);
    
    %% Part C 
    
    % This approach is memory (space) efficient as it resues the same scalars
    % and doesn't store them anywhere, so it will never grow in size. Thus
    % it has a memory efficiency of O(1). However, since this is an iterative
    % process by the number of trials, it will have O(N) time complexity.
    
    N = 1e5;
    map_errors = 0;
    ml_errors = 0;
    
        for i = 1:N
            [u1_color, u2_color] = select_2(r1, b1, r2, b2);
    
            map_decision = map_decision_rule(u2_color);
            ml_decision  = ml_decision_rule(u2_color);
            
            % incrememnt error by 1 if result is not the same as predicted
            if map_decision ~= u1_color
                map_errors = map_errors + 1;
            end
            if ml_decision ~= u1_color
                ml_errors = ml_errors + 1;
            end
        end
    
        ex_map_error = map_errors / N;
        ex_ml_error = ml_errors / N;
    
        fprintf('\nTheoretical MAP Error Probability: %f\n', (map_error_Red2 + map_error_Blue2));
        fprintf('Experimental MAP Error Probability: %f\n', ex_map_error);
        fprintf('\nTheoretical ML Error Probability: %f\n', (ml_error_Red2 + ml_error_Blue2));
        fprintf('Experimental ML Error Probability: %f\n', ex_ml_error);
end


function decision_rule = MAP_decision(Posterior)
    decision_rule = zeros(1, 2);
    
    % Red2
    % Choose largest posterior P(Red1 | Red2) or P(Blue1 | Red2)
    if Posterior(1, 1) > Posterior(1, 2)
        decision_rule(1) = 1; % Red1
    else
        decision_rule(1) = 2; % Blue1
    end
    
    % Blue2
    % Choose largest posterior P(Red1 | Blue2) or P(Blue1 | Blue2)
    if Posterior(2, 1) > Posterior(2, 2)
        decision_rule(2) = 1; % Red1
    else
        decision_rule(2) = 2; % Blue1
    end

    fprintf('\nMAP Decision Rule: [%d, %d]\n\n', decision_rule);
end


function decision_rule = ML_decision(likelihoods)
    decision_rule = zeros(1, 2);

    % Red2
    % Choose largest likelihoods P(Red2 | Red1) or P(Red2 | Blue1)
    if likelihoods(1, 1) > likelihoods(2, 1)
        decision_rule(1) = 1; % Red1
    else
        decision_rule(1) = 2; % Blue1
    end

    % Blue2
    % Choose largest likelihood P(Blue2 | Red1) or P(Blue2 | Blue1)
    if likelihoods(1, 2) > likelihoods(2, 2)
        decision_rule(2) = 1; % Red1
    else
        decision_rule(2) = 2; % Blue1
    end

    fprintf('\nML Decision Rule: [%d, %d]\n\n', decision_rule);
end


% joint error probabilities for each outcome (Red2 and Blue2)
function [err_r2, err_b2] = error_prob(decision_rule, P_prior, likelihoods)
    
% Red2 observed
    if decision_rule(1) == 1 
        % Decision says Red1 on Red2. Wrong if true class was Blue1.
        % P(Red2 | Blue1) * P(Blue1)
        err_r2 = likelihoods(2,1) * P_prior(2); 
    else 
        % Decision says Blue1 on Red2. Wrong if true class was Red1.
        % P(Red2 | Red1) * P(Red1)
        err_r2 = likelihoods(1,1) * P_prior(1); 
    end

    % Blue2 observed
    if decision_rule(2) == 1 
        % Decision says Red1 on Blue2. Wrong if true class was Blue1
        % P(Blue2 | Blue1) * P(Blue1)
        err_b2 = likelihoods(2,2) * P_prior(2); 
    else 
        % Decision says Blue1 on Blue2. Wrong if true class was Red1.
        % P(Blue2 | Red1) * P(Red1)
        err_b2 = likelihoods(1,2) * P_prior(1); 
    end
end


%% Part 1 - select from 1 urn
function urn1color = select_1(r1, b1)
    total = r1 + b1;
    picked = randi(total);

    if picked <= r1
        urn1color = 1; % Red 
    else 
        urn1color = 2; % Blue
    end
end


%% Part 2 - select from one urn, then add to the first ball to the second urn and reselect
function [u1_color, u2_color] = select_2(r1, b1, r2, b2)
    u1_total = r1 + b1;
    pick1 = randi(u1_total);
    
    if pick1 <= r1
        u1_color = 1; % Red 
        r2 = r2 + 1; 
    else
        u1_color = 2; % Blue
        b2 = b2 + 1;
    end

    u2_total = r2 + b2;
    pick2 = randi(u2_total);

    if pick2 <= r2
        u2_color = 1; % Red 
    else
        u2_color = 2; % Blue
    end
end

