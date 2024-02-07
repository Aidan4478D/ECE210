% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment03.m -- "Amogus"
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 

ITERATIONS = 1e1;

CREWMATES = 6;
ROUNDS = 12;

CREWMATE_SIDES = 4;
IMPOSTER_ROLLS = 2;
IMPOSTER_SIDES = 2;

rng(0x73757300); % random number seeding function

%% part 1

% randi(bound, dim, iterations)

% each crewmate's sus resistance
% one roll of a four sided die
crewmates = randi(CREWMATE_SIDES, CREWMATES, ITERATIONS);

% imposter's sus ability
% sum of two rolls of a two sided die
sus = sum(randi(IMPOSTER_SIDES, IMPOSTER_ROLLS, ITERATIONS));

% twelve rolls of six sided die will be the targets for each game
targets = randi(CREWMATES, ROUNDS, ITERATIONS);


%% part 2

% kills matrix represents crewmates that are alive and dead

% start with each crewmate being alive (ordered 1-6 by column with values 1 = alive and 0 = dead)
% initial kills array should be array of 1's as no resistance is <= 0
kills = logical(crewmates);

% can_die = zeros(size(crewmates));
% died = crewmates(1:CREWMATES, :) < sus

% col = targets(1:ROUNDS, :)


