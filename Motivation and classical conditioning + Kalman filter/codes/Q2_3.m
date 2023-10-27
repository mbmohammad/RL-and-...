%% set params
clc;clear all;close all;
%%
n_trials = 60;
r = 1;

InitialVar = 0.6;
measurementNoiseVar = 0.1;
ProcessNoiseVar = 0.004;

Weight = zeros(2, n_trials);
sigma = zeros(2, 2, n_trials);
sigma_predict = zeros(2, 2, n_trials);

W = [ProcessNoiseVar 0; 0 ProcessNoiseVar];
V = measurementNoiseVar;

sigma(:,:, 1) = [InitialVar 0; 0 InitialVar];
Weight(:, 1) = [0 0];

input = [ones(1, n_trials); zeros(1, n_trials)];
reward = r*ones(1, n_trials);
reward(n_trials/2+1:end) = -reward(n_trials/2+1:end);
%% LOOP
for i = 1 : n_trials-1
    sigma_predict(:,:, i+1) = sigma(:, :, i) + W;
%     sigma_predict(:,:, i+1)*input(:, i)'
%     input(:, i)*sigma_predict(:,:, i+1)*input(:, i)'
    G = sigma_predict(:,:, i+1)*input(:, i)/(input(:, i)'*sigma_predict(:,:, i+1)*input(:, i)+V);
    sigma(:,:, i+1) = sigma_predict(:,:, i+1) - G*input(:, i)'*sigma_predict(:,:, i+1);
    Weight(:, i+1) = Weight(:, i) + G*(reward(:, i) - input(:, i)'*Weight(:, i));
    if(i == n_trials/4)
        sigma(2,2, i+1) = InitialVar;
    end
end
%%
figure
plot(Weight(1,:))
xlabel("Trial", 'interpreter', 'Latex')
ylabel("W", 'interpreter', 'Latex')
title(["Paradigm" , sprintf("Process Noise Variance = %.3f", ProcessNoiseVar), sprintf("Measurement Noise Variance = %.3f ", measurementNoiseVar)])

figure
plot(reshape(sigma(1,1, :),1,n_trials))
xlabel("Trial", 'interpreter', 'Latex')
ylabel("Variance", 'interpreter', 'Latex')
title(["Paradigm" , sprintf("Process Noise Variance = %.3f", ProcessNoiseVar), sprintf("Measurement Noise Variance = %.3f ", measurementNoiseVar)])


