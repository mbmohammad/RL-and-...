%% set parameters
clc;clear all;close all;
epsilon = 0.004;
r = 1;
n_trials = 2000;
S1 = ones(1, n_trials);
S2 = [zeros(1, n_trials/2) ones(1, n_trials/2)];
S = [S1;S2];
R1 = r*ones(1, n_trials);
% W = 0.5*ones(2, n_trials);
W = zeros(2, n_trials);

%% process
for i = 2 : n_trials
    W(:, i) = W(:, i-1)+epsilon*S(:, i-1)*(R1(i-1)-W(1, i-1)*S(1,i-1) - W(2, i-1)*S(2,i-1));
end
%% plot
figure
plot(W(1,:), 'color', 'r', 'LineWidth', 2)
hold on
plot(W(2,:), 'color', 'g', 'LineWidth', 2)
xlabel("Trial", 'interpreter', 'Latex')
ylabel("W", 'interpreter', 'Latex')
title(sprintf("Blocking Paradigm, r = %.2f, epsilon = %.3f", r, epsilon), 'interpreter', 'Latex')
legend("W1", "W2")



