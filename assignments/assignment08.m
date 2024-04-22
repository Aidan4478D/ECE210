% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment08.m -- Under Pressure
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 

%% part 1
syms Ic Vbe Vt Ae q Dn ni Nb Wb;

Is = Ae * q * Dn * ni^2 / (Nb * Wb);
eq = Ic == Is * exp(Vbe / Vt);

solve(eq, Ae)


%% part 2
syms Vs(t) R C Vc(t);

eq = Vs(t) == R * C * diff(Vc, t) + Vc(t);
dsolve(eq)

%% part 3

mu0 = sym(4 * 1e-7);
mu0_f = vpa(mu0 * pi, digits(floor(1000 * pi)))
