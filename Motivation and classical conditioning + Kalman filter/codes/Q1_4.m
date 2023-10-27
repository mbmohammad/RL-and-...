%% set parameters
clc;clear all;close all;
epsilon = 0.03;
r = 1;
alpha = 0.5;
n_trials = 600;
S1 = ones(1, n_trials);
S2 = rand(1, n_trials);
S2 = double(S2<alpha);
S = [S1;S2];
R1 = r*ones(1, n_trials);
R1(S2==1) = 0;
W = zeros(2, n_trials);
W(:,1) = 2*[0.85 0.85];
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
title(sprintf("Inhibitory Paradigm, r = %.2f, epsilon = %.2f", r, epsilon), 'interpreter', 'Latex')
legend("W1", "W2")


