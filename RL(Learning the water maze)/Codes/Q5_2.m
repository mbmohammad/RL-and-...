clc; clear all;close all;
% Initialize the map (15 x 15 grid)
map = zeros(15, 15);

% Define the goal state
goal_state1 = [5,5];
goal_state2 = [10,10];

% Define the reward matrix (initialize to -1)
reward_matrix = -1 * ones(15, 15);
reward_matrix(5, 5) = 10; % When reaching the goal state, obtain a reward of 10
reward_matrix(10, 10) = -10;
% Define the Q-matrix
qmatrix = zeros(15, 15, 4);

% Define the hyperparameters
num_episodes = 10000;
gamma = 0.8; % discount factor
epsilon = 0.1; % epsilon-greedy exploration rate
lambda = 0.5;

heatMap = zeros(21,20);

hm1 = 0;
hm2 = 0;
for  alpha = 0.2:0.05:0.5
    hm1  = hm1 + 1
    hm2=0;
    for gamma = 0.5:0.05:1
        gamma
        hm2  = hm2 + 1;
        % Define the Q-matrix
qmatrix = zeros(15, 15, 4);
numberOfActions = zeros(1, num_episodes);
% Q-learning algorithm
for i = 1:num_episodes
    % Define the starting state randomly
    state = [randi(15), randi(15)];
    k = 0;%number of moves
    e = zeros(15,15);
    % Loop until reaching the goal state
    while ~(isequal(state, goal_state1) || isequal(state, goal_state2))
        k = k+1;
        % Choose an action using epsilon-greedy exploration strategy
        if rand() < epsilon
            action = randi(4);
        else
            [~, action] = max(qmatrix(state(1), state(2), :));
        end
        
        % Compute the next state and reward
        next_state = getNextState(state, action);
        reward = reward_matrix(next_state(1), next_state(2));
        e(state(1), state(2)) = e(state(1), state(2))+1;
        % Update the Q-matrix using the Q-learning update rule
        qmatrix(state(1), state(2), action) = qmatrix(state(1), state(2), action)...
            + e(state(1), state(2))*alpha *(reward + gamma * max(qmatrix(next_state(1), next_state(2), :))- qmatrix(state(1), state(2), action));
       
        e(state(1), state(2)) = e(state(1), state(2))* gamma * lambda;
        % Move to the next state
        state = next_state;
    end
    numberOfActions(i)=k;
end
    for s = 40 : num_episodes
        if(mean(numberOfActions(s-39:s))<14)
            heatMap(hm1,hm2) = s;
            break;
        end
    end
    end
end
%%
figure
heatmap(0.5:0.05:1,0.2:0.05:0.5, heatMap(1:7,1:11))
title("minimum number of iterations needed for reaching the optimal paths heat map, lambda = 0.5")
xlabel("discount factor")
ylabel("learning rate")
%%
figure
plot(0.2:0.05:0.5, mean(heatMap(1:7,1:11),2),'LineWidth',2,'Color','r')
title(["minimum number of iterations needed for reaching" "the optimal paths heat map, $\lambda = 0.5$(mean over discount factors)"], 'interpreter','latex')
xlabel("learning rate",'interpreter','latex')
ylabel("number of iteration",'interpreter','latex')
%%
figure
plot(0.5:0.05:1, mean(heatMap(1:7,1:11),1),'LineWidth',2,'Color','r')
title(["minimum number of iterations needed for reaching" "the optimal paths heat map, $\lambda = 0.5$(mean over learning rates)"], 'interpreter','latex')
xlabel("discount factor",'interpreter','latex')
ylabel("number of iteration",'interpreter','latex')
%% Function to compute the next state given the current state and action
function next_state = getNextState(curr_state, action)
    switch action
        case 1 % Move up
            next_state = [max(curr_state(1)-1, 1), curr_state(2)];
        case 2 % Move right
            next_state = [curr_state(1), min(curr_state(2)+1, 15)];
        case 3 % Move down
            next_state = [min(curr_state(1)+1, 15), curr_state(2)];
        case 4 % Move left
            next_state = [curr_state(1), max(curr_state(2)-1, 1)];
    end
end
