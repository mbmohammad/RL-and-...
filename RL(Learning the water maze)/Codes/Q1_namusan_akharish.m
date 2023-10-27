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
alpha = 0.5; % learning rate
gamma = 0.8; % discount factor
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
%% prepare selected episodes and image them
pauses = [0.1 0.1 0.2*ones(1, length(recorded_episodes)-2)];
for l = 1 : length(recorded_episodes)
selected_episodes = zeros(2, 1000);
selected_index = l;
for i = 1 : 1000
   selected_episodes(:,i) = episode_recorder(:, selected_index,i); 
end
selected_episodes( :, ~any(selected_episodes,1) ) = [];
if(sum(abs(selected_episodes(:,end)-[5;5]))<=2)
    selected_episodes = [selected_episodes [goal_state1(1);goal_state1(2)]];
else
    selected_episodes = [selected_episodes [goal_state2(1);goal_state2(2)]];
end

figure
    for i = 1 : size(selected_episodes,2)
    map = 100*ones(15,15);
    map(goal_state1(1),goal_state1(2)) = 140;
    map(goal_state2(1),goal_state2(2)) = 220;
    map(selected_episodes(1,i),selected_episodes(2,i)) = 165;
    image(map);colormap(jet(256));
    title(sprintf("Path of rat in trial %d", recorded_episodes(selected_index)));
    xlabel("X")
    ylabel("Y")
    pause(pauses(l))
    hold on
    end
end
%% plot path
for l = 1 : length(recorded_episodes)
    selected_episodes = zeros(2, 1000);
selected_index = l;
for i = 1 : 1000
   selected_episodes(:,i) = episode_recorder(:, selected_index,i); 
end
selected_episodes( :, ~any(selected_episodes,1) ) = [];
if(sum(abs(selected_episodes(:,end)-[5;5]))<=2)
    selected_episodes = [selected_episodes [goal_state1(1);goal_state1(2)]];
else
    selected_episodes = [selected_episodes [goal_state2(1);goal_state2(2)]];
end
    figure;
    plot(selected_episodes(1,:), selected_episodes(2,:), '--b','LineWidth',2)
    xlim([1 15])
    ylim([1 15])
    xlabel("X")
    ylabel("Y")
    axis ij
    hold on
    scatter(5,5,'g','filled')
    hold on
    scatter(10,10,'r','filled')
    hold on
    scatter(selected_episodes(1,1), selected_episodes(2,1),'y','filled')
    title(sprintf("Rat Path in trial %d", recorded_episodes(selected_index)));
    legend('rat path','reward','cat','starting point')
end


%%
figure
V_Y = -qmatrix(:,:, 1) + qmatrix(:,:, 3);
V_X = qmatrix(:,:, 2) - qmatrix(:,:, 4);
quiver(V_X,V_Y,'r')
hold on
contour(sqrt(V_X.^2+V_Y.^2))
title("Gradients and Contour plot of the learned values gamma = 0.85 alpha = 0.6", 'interpreter','latex')
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
