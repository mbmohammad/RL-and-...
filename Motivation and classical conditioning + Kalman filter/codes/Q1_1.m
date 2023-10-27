%% set parameters
clc;clear all;close all;
r = 0.8;
epsilon = 0.1;
n_trials = 200;
S1 = ones(1, n_trials);
R1 = r*[ones(1, n_trials/2) zeros(1, n_trials/2)];
W = zeros(1, n_trials);
%% process
for i = 2 : n_trials
    W(i) = W(i-1)+epsilon*S1(i-1)*(R1(i-1)-W(i-1)*S1(i-1));
end
%% plot
figure
plot(W, 'color', 'r', 'LineWidth', 2)
xlabel("Trial", 'interpreter', 'Latex')
ylabel("W", 'interpreter', 'Latex')
title(sprintf("Extinction Paradigm, r = %.2f, epsilon = %.2f", r, epsilon), 'interpreter', 'Latex')


