%% set parameters
clc;clear all;close all;
epsilon = 0.01;
alpha = 0.3;
n_trials = 1000;
S1 = ones(1, n_trials);
R1 = rand(n_trials, 1);
R1 = double(R1<alpha);
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
title("Partial Paradigm", 'interpreter', 'Latex')
title(sprintf("Partial Paradigm, alpha(density of ones) = %.2f, epsilon = %.2f", alpha, epsilon), 'interpreter', 'Latex')

