% SPDX-License-Identifier: GPL-3.0-or-later
%
% ECE210 assignment05.m -- Plotting, Scheming Even
% Copyright (C) 2024 Aidan Cusa <aidancusa@gmail.com>

clc;    % clear command window
clear;  % clear all variables from current worwkspace
close all; 


%% part 1
n = 50;
t = linspace(-pi, pi, 1000);

a_n = (2 .* (0:n) + 1)';

terms = zeros(length(t), n); % matrix for individual terms
sq = sum(sin(a_n .* t) ./ a_n); % THE square wave

for i = 1:n
    terms(:, i) = sin(a_n(i) * t) ./ a_n(i); % fill matrix with terms
end
                   
figure; 
hold on;

plot(t, terms)
plot(t, sq)

title('Square Wave with n=50');
xlabel('X [rad]');
xlim([t(1), t(end)])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
xticks([-pi, -pi/2, 0, pi/2, pi]);
ylabel('y');

%% part 2
figure;
subplot(2, 1, 1);
plot(t, sq)
%plot(t, fourier_sum(0:50, t))
title("Square Wave n=50 Approximation");
xlabel('X [rad]');
xlim([t(1), t(end)])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
xticks([-pi, -pi/2, 0, pi/2, pi]);
ylabel('y');

subplot(2, 1, 2);
plot(t, terms)
title("Square Wave n= 0-49 Approximation");
xlabel('X [rad]');
xlim([t(1), t(end)])
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
xticks([-pi, -pi/2, 0, pi/2, pi]);
ylabel('y');

hold on;

sgtitle("Fourier Square Wave Approximation");
hold off; % is this necessary?

%% part 3

S = linspace(-2*pi, 2*pi);
x = S;
y = S; 

[X, Y] = meshgrid(x, y);
z = X .* sin(X) - Y .* cos(Y);

figure;
surf(X, Y, z)

%% extra matlab logo

% lowkey cheated with the link below but it still looks cool 
% https://www.mathworks.com/help/matlab/visualize/creating-the-matlab-logo.html

L = 160*membrane(1,100);
figure;
logo = surf(L);
logo.EdgeColor = 'none';

l1 = light;
l1.Position = [160 400 80];
l1.Style = 'local';
l1.Color = [0 0.8 0.8];
 
l2 = light;
l2.Position = [.5 -1 .4];
l2.Color = [0.8 0.8 0];

logo.FaceColor = [0.9 0.2 0.2];
logo.FaceLighting = 'gouraud';
logo.AmbientStrength = 0.3;
logo.DiffuseStrength = 0.6; 
logo.BackFaceLighting = 'lit';

logo.SpecularStrength = 1;
logo.SpecularColorReflectance = 1;
logo.SpecularExponent = 7;

axis off
