%% set params
clc;clear all;close all;
%%
Gs = zeros(20, 20, 2);
mNoise = 0.05:0.05:1;
pNoise = 0.005:0.005:0.1;
X1 = zeros(20, 20);
X2 = zeros(20, 20);


for k1 = 1 : 20
    for k2 = 1 : 20
n_trials = 36;
r = 1;

InitialVar = 0.6;
measurementNoiseVar = mNoise(k1);
ProcessNoiseVar = pNoise(k2);

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
    Gs(k1,k2,:) = G;
    X1(k1,k2) = mNoise(k1);
    X2(k1,k2) = pNoise(k1);
    end
end

%%
figure
for i = 1 : 20
plot(mNoise, Gs(:,i,1))
xlabel("measurement noise", 'interpret', 'latex')
title(["Kalman Filter Gain for first stimulus", "with different levels of process noise"], 'interpret', 'latex')
ylabel("Gain", 'interpret', 'latex')
hold on
end


figure
for i = 1 : 20
plot(pNoise, Gs(i,:,1))
xlabel("process noise", 'interpret', 'latex')
title(["Kalman Filter Gain for first stimulus", "with different levels of process noise"], 'interpret', 'latex')
ylabel("Gain", 'interpret', 'latex')
hold on
end





