% %% clear 
% clc;clear all;close all;
% %% define W and R
% n_trials = 100;
% W = zeros(1, n_trials);
% V_mu = 0;
% V_sigma = 0.01;
% C = zeros(1,n_trials);
% C(40) = 1;C(80) = 1;
% PHI_mu = 0;
% PHI_sigma = 4;
% 
% for i = 2 : length(W)
%     W(i) = W(i-1)+ normrnd(V_mu,V_sigma) + C(i)*normrnd(PHI_mu,PHI_sigma);
% end
% 
% R = zeros(1, n_trials);
% X = ones(1, n_trials);
% Eta_mu = 0;
% Eta_sigma = 0.35;
% 
% for i = 1 : length(R)
%     R(i) = W(i) + normrnd(Eta_mu,Eta_sigma);
% end
%%
figure
plot(W, 'b')
hold on
plot(R, 'r*')
legend('W(t)', 'R(t)', 'interpreter', 'latex', 'location', 'northEast')
xlabel("Trials", 'interpreter', 'latex')
title("Sample sequence of w(t) and corresponding sequence of noisy r(t)", 'interpreter', 'latex')

%% find W_hat
InitialVar = 0.6;
measurementNoiseVar = Eta_sigma;
ProcessNoiseVar = V_sigma;

Weight = zeros(2, n_trials);
sigma = zeros(2, 2, n_trials);
sigma_predict = zeros(2, 2, n_trials);

W1 = [ProcessNoiseVar 0; 0 ProcessNoiseVar];
V = measurementNoiseVar;

sigma(:,:, 1) = [InitialVar 0; 0 InitialVar];
Weight(:, 1) = [0 0];

reward=R;
input = [X; zeros(1, n_trials)];

beta = zeros(1, n_trials);
Thr = 3;
for i = 1 : n_trials-1
    sigma_predict(:,:, i+1) = sigma(:, :, i) + W1;
%     sigma_predict(:,:, i+1)*input(:, i)'
%     input(:, i)*sigma_predict(:,:, i+1)*input(:, i)'
    G = sigma_predict(:,:, i+1)*input(:, i)/(input(:, i)'*sigma_predict(:,:, i+1)*input(:, i)+V);
    sigma(:,:, i+1) = sigma_predict(:,:, i+1) - G*input(:, i)'*sigma_predict(:,:, i+1);
    Weight(:, i+1) = Weight(:, i) + G*(reward(:, i) - input(:, i)'*Weight(:, i));
    if(i == n_trials/4)
        sigma(2,2, i+1) = InitialVar;
    end
    beta(i) = (reward(:, i) - input(:, i)'*Weight(:, i))^2/(input(:, i)'*sigma_predict(:,:, i+1)*input(:, i)+V);
%     if(i>1)
%     if(beta(i)>Thr)
%         sigma(:,:, i+1) = sigma(:,:, i+1) + 50;
%     end
%     end
end


figure
plot(W, 'b')
hold on
plot(R, 'r*')
hold on
plot(Weight(1, :),'co')
legend('W(t)', 'R(t)','$\hat{W}$', 'interpreter', 'latex', 'location', 'northWest')
xlabel("Trials", 'interpreter', 'latex')
title("Sample sequence of w(t) and corresponding sequence of noisy r(t) and predicted W by ordinary Kalman filter", 'interpreter', 'latex')

%% find W_hat
InitialVar = 0.6;
measurementNoiseVar = Eta_sigma;
ProcessNoiseVar = V_sigma;

Weight = zeros(2, n_trials);
sigma = zeros(2, 2, n_trials);
sigma_predict = zeros(2, 2, n_trials);

W1 = [ProcessNoiseVar 0; 0 ProcessNoiseVar];
V = measurementNoiseVar;

sigma(:,:, 1) = [InitialVar 0; 0 InitialVar];
Weight(:, 1) = [0 0];

reward=R;
input = [X; zeros(1, n_trials)];

beta = zeros(1, n_trials);
Thr = 2;
for i = 1 : n_trials-1
    sigma_predict(:,:, i+1) = sigma(:, :, i) + W1;
%     sigma_predict(:,:, i+1)*input(:, i)'
%     input(:, i)*sigma_predict(:,:, i+1)*input(:, i)'
    G = sigma_predict(:,:, i+1)*input(:, i)/(input(:, i)'*sigma_predict(:,:, i+1)*input(:, i)+V);
    sigma(:,:, i+1) = sigma_predict(:,:, i+1) - G*input(:, i)'*sigma_predict(:,:, i+1);
    Weight(:, i+1) = Weight(:, i) + G*(reward(:, i) - input(:, i)'*Weight(:, i));
    if(i == n_trials/4)
        sigma(2,2, i+1) = InitialVar;
    end
    beta(i) = (reward(:, i) - input(:, i)'*Weight(:, i))^2/(input(:, i)'*sigma_predict(:,:, i+1)*input(:, i)+V);
    if(i>1)
    if(beta(i)>Thr)
        sigma(:,:, i+1) = sigma(:,:, i+1) + 50;
    end
    end
end


figure
plot(W, 'b')
hold on
plot(R, 'r*')
hold on
plot(Weight(1, :),'co')
legend('W(t)', 'R(t)','$\hat{W}$', 'interpreter', 'latex', 'location', 'northWest')
xlabel("Trials", 'interpreter', 'latex')
title("Sample sequence of w(t) and corresponding sequence of noisy r(t) and predicted W by modified Kalman filter", 'interpreter', 'latex')

%%

ACh = 12*double(beta>Thr);

figure
plot(beta)
hold on
yline(Thr, '--r')
hold on
plot(ACh, '-.m')
legend('NE', '$\gamma$','ACh', 'interpreter', 'latex', 'location', 'northEast')
xlabel("Trials", 'interpreter', 'latex')
title("Corresponding levels of ACh (dashed) and NE (solid) for the same sequence", 'interpreter', 'latex')


sum((R-Weight(1, :)).^2)










