%% set parameters
clc;clear all;close all;
epsilon = 0.001;
alpha1 = 0.4;
alpha2 = 0.2;
r = 1;
n_trials = 5000;
S1 = rand(n_trials, 1);
S1 = double(S1<alpha1)';
S2 = rand(n_trials, 1);
S2 = double(S2<alpha2)';
S = [S1;S2];
R1 = r*ones(1, n_trials);
R1(S1==0 & S2==0) = 0;
R1((S1==1 & S2==0) | (S1==0 & S2==1)) = 0;
W = zeros(2, n_trials);
%% process
for i = 2 : length(R1)
    W(:, i) = W(:, i-1)+epsilon*S(:, i-1)*(R1(i-1)-W(1, i-1)*S(1,i-1) - W(2, i-1)*S(2,i-1));
end
%% plot
figure
plot(W(1,:), 'color', 'r', 'LineWidth', 2)
hold on
plot(W(2,:), 'color', 'g', 'LineWidth', 2)
xlabel("Trial", 'interpreter', 'Latex')
ylabel("W", 'interpreter', 'Latex')
title(sprintf("Overshadow Paradigm, r = %.2f, alpha1 = %.2f, alpha2 = %.2f, epsilon= %.3f", r, alpha1, alpha2, epsilon), 'interpreter', 'Latex')
legend("W1", "W2")

