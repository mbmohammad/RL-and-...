%% set params
clc;clear all;close all;
%%
n_trials = 60;
r = 1;

InitialVar = 0.6;
measurementNoiseVar = 0.1;
ProcessNoiseVar = 0.002;

Weight = zeros(2, n_trials);
sigma = zeros(2, 2, n_trials);
sigma_predict = zeros(2, 2, n_trials);

W = [ProcessNoiseVar 0; 0 ProcessNoiseVar];
V = measurementNoiseVar;

sigma(:,:, 1) = [InitialVar 0; 0 InitialVar];
Weight(:, 1) = [0 0];

input = ones(2, n_trials);
input(2, 1:n_trials/4) = 0; %blocking
reward = r*ones(1, n_trials);
%reward(n_trials/4+1:n_trials) = 2*reward(n_trials/4+1:n_trials); %for unblocking
% input(2, n_trials/4+1:end) = 0; %backward blocking
%% part e
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
hold on
plot(n_trials/4+1:n_trials, Weight(2,n_trials/4+1:end))
xlabel("Trial", 'interpreter', 'Latex')
ylabel("W", 'interpreter', 'Latex')
title(["Backward Blocking Paradigm" , sprintf("Process Noise Variance = %.3f", ProcessNoiseVar), sprintf("Measurement Noise Variance = %.3f ", measurementNoiseVar)])
xline(n_trials/4+1, '--g',{'S2 is presented from here'})
legend("$\hat{W_{1}}(t)$","$\hat{W_{2}}(t)$", 'interpreter', 'latex')

figure
plot(reshape(sigma(1,1, :),1,n_trials))
hold on
plot(n_trials/4+1:n_trials, reshape(sigma(2,2,n_trials/4+1:end),1, floor(3*n_trials/4)))
xlabel("Trial", 'interpreter', 'Latex')
ylabel("Variance", 'interpreter', 'Latex')
xline(n_trials/4+1, '--g',{'S2 is presented from here'})
title(["Backward Blocking Paradigm" , sprintf("Process Noise Variance = %.3f", ProcessNoiseVar), sprintf("Measurement Noise Variance = %.3f ", measurementNoiseVar)])
legend("$\sigma_{1}^2(t)$","$\sigma_{2}^2(t)$", 'interpreter', 'latex')


%% contour plots
m = [Weight(1,1) Weight(2,1)];
S = sigma(:,:,1);
figure
G = gmdistribution(m,S);
F = @(x,y) pdf(G,[x y]);
fcontour(F, [-2 2 -2 2])
xlabel("$\bar{W_1}$", 'interpreter', 'Latex')
ylabel("$\bar{W_2}$", 'interpreter', 'Latex')
title(["Backward Blocking Paradigm trial 1" , sprintf("Process Noise Variance = %.3f", ProcessNoiseVar), sprintf("Measurement Noise Variance = %.3f ", measurementNoiseVar)])

m = [Weight(1,9) Weight(2,9)];
S = sigma(:,:,9);
figure
G = gmdistribution(m,S);
F = @(x,y) pdf(G,[x y]);
fcontour(F, [-2 2 -2 2])
xlabel("$\bar{W_1}$", 'interpreter', 'Latex')
ylabel("$\bar{W_2}$", 'interpreter', 'Latex')
title(["Backward Blocking Paradigm trial 9" , sprintf("Process Noise Variance = %.3f", ProcessNoiseVar), sprintf("Measurement Noise Variance = %.3f ", measurementNoiseVar)])

m = [Weight(1,19) Weight(2,19)];
S = sigma(:,:,19);
figure
G = gmdistribution(m,S);
F = @(x,y) pdf(G,[x y]);
fcontour(F, [-2 2 -2 2])
xlabel("$\bar{W_1}$", 'interpreter', 'Latex')
ylabel("$\bar{W_2}$", 'interpreter', 'Latex')
title(["Backward Blocking Paradigm trial 19" , sprintf("Process Noise Variance = %.3f", ProcessNoiseVar), sprintf("Measurement Noise Variance = %.3f ", measurementNoiseVar)])





