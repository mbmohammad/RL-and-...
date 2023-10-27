clc; clear all;close all;
% Initialize the map (15 x 15 grid)
map = zeros(15, 15);

% Define the goal state
goal_state1 = [5,5];
goal_state2 = [10,10];
goal_state3 = [13,13];


% Define the reward matrix (initialize to -1)
reward_matrix = -1 * ones(15, 15);
reward_matrix(5, 5) = 10; % When reaching the goal state, obtain a reward of 10
reward_matrix(10, 10) = -10;
reward_matrix(13, 13) = 30;
% Define the Q-matrix
qmatrix = zeros(15, 15, 4);

% Define the hyperparameters
num_episodes = 10000;
alpha = 0.1; % learning rate
gamma = 0.65; % discount factor
epsilon = 0.1; % epsilon-greedy exploration rate

numberOfActions = zeros(1, num_episodes);
episode_recorder = zeros(2, 11, 1000);
recorded_episodes = [10, 30, 60, 120, 240, 700,996,997,998, 999, 1000];
%% Q-learning algorithm
for i = 1:num_episodes
    % Define the starting state randomly
    state = [randi(15), randi(15)];
    k = 0;%number of moves
    % Loop until reaching the goal state
    while ~(isequal(state, goal_state1) || isequal(state, goal_state2) || isequal(state, goal_state3))
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
        
        % Update the Q-matrix using the Q-learning update rule
        qmatrix(state(1), state(2), action) = qmatrix(state(1), state(2), action)...
            + alpha *(reward + gamma * max(qmatrix(next_state(1), next_state(2), :))- qmatrix(state(1), state(2), action));
        
        for l = 1 : length(recorded_episodes)
            if(i == recorded_episodes(l))
                episode_recorder(:, l, k) = state;
            end
        end
        
        % Move to the next state
        state = next_state;
    end
    numberOfActions(i)=k;
end
%%
figure
V_Y = -qmatrix(:,:, 1) + qmatrix(:,:, 3);
V_X = qmatrix(:,:, 2) - qmatrix(:,:, 4);
quiver(V_X,V_Y,'r')
hold on
contour(sqrt(V_X.^2+V_Y.^2))
title("Gradients and Contour plot of the learned values gamma = 0.65 alpha = 0.1", 'interpreter','latex')
xlabel("X", 'interpreter','latex')
ylabel("Y", 'interpreter','latex')
axis ij

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
