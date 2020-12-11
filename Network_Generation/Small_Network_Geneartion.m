%Random Small Network Generation
%Author: Toru kitagawa and Guanyi Wang

clear all
filename = '../simulation_small.mat';
%generate random graph
%randomly generate adjacency matrix
rng(0)%set seed = 0
% for N =25
N_1 = 25; %first choice of nodes, N = 25

e_1 = 7*N_1; %first choice of edges = 7N 
adj_1 = random_graph(N_1,[],e_1);%generate random adjacency matrix given number of nodes and edges

e_2 = 10*N_1; %second choice of edges = 10N
adj_2 = random_graph(N_1,[],e_2);

e_3 = fulledges(N_1);  %third choice of edges = full
adj_3 = random_graph(N_1,[],e_3);

%generate random state variable for N =25
P_s = 0.4;%probability of being susceptible at initial stage
P_i = 0.5;%probability of being infected at initial stage
idx_state_1 = randperm(N_1);
S_1 = zeros(N_1,1);
S_1(idx_state_1(1:round(P_s*N_1)),:) = 1;%randomly generate the susceptible units with p_s
I_1 = zeros(N_1,1);
I_1(idx_state_1(round(P_s*N_1)+1:round(P_s*N_1+P_i*N_1)),:) = 1;%randomly generate the infected units with p_i
R_1 = zeros(N_1,1);
R_1(idx_state_1(round(P_s*N_1+P_i*N_1)+1:end),:) = 1;%randomly generate the recovered units with p_s

%random allocate group for N = 25
P_y = 0.4;%probability of staying group 1
idx_group_1 = randperm(N_1);
a_1 = zeros(N_1,1);
a_1(idx_group_1(1:round(P_y*N_1)),:) = 1;%randomly generate the group 1 units with p_y
b_1 = zeros(N_1,1);
b_1(idx_group_1(round(P_y*N_1)+1:end),:) = 1;%randomly generate the group 2 units with (1-p_y)

% for N =30
N_2 = 30;%second choice of nodes, N = 30
e_4 = 7*N_2; 
adj_4 = random_graph(N_2,[],e_4);
e_5 = 10*N_2; 
adj_5 = random_graph(N_2,[],e_5);
e_6 = fulledges(N_2);
adj_6 = random_graph(N_2,[],e_6);



%generate random state variable for N = 30
idx_state_2 = randperm(N_2);
S_2 = zeros(N_2,1);
S_2(idx_state_2(1:round(P_s*N_2)),:) = 1;
I_2 = zeros(N_2,1);
I_2(idx_state_2(round(P_s*N_2)+1:round(P_s*N_2+P_i*N_2)),:) = 1;
R_2 = zeros(N_2,1);
R_2(idx_state_2(round(P_s*N_2+P_i*N_2)+1:end),:) = 1;

%random allocate group for N = 30
idx_group_2 = randperm(N_2);
a_2 = zeros(N_2,1);
a_2(idx_group_2(1:round(P_y*N_2)),:) = 1;
b_2 = zeros(N_2,1);
b_2(idx_group_2(round(P_y*N_2)+1:end),:) = 1;

% for N =35
N_3 = 35;%third choice of nodes, N = 35
e_7 = 7*N_3; 
adj_7 = random_graph(N_3,[],e_7);
e_8 = 10*N_3; 
adj_8 = random_graph(N_3,[],e_8);
e_9 = fulledges(N_3);
adj_9 = random_graph(N_3,[],e_9);

%generate random state variable for N = 35
idx_state_3 = randperm(N_3);
S_3 = zeros(N_3,1);
S_3(idx_state_3(1:round(P_s*N_3)),:) = 1;
I_3 = zeros(N_3,1);
I_3(idx_state_3(round(P_s*N_3)+1:round(P_s*N_3+P_i*N_3)),:) = 1;
R_3 = zeros(N_3,1);
R_3(idx_state_3(round(P_s*N_3+P_i*N_3)+1:end),:) = 1;

%random allocate group for N = 35
idx_group_3 = randperm(N_3);
a_3 = zeros(N_3,1);
a_3(idx_group_3(1:round(P_y*N_3)),:) = 1;
b_3 = zeros(N_3,1);
b_3(idx_group_3(round(P_y*N_3)+1:end),:) = 1;

% Save coefficients and other variables for graphs
save(filename,'S_1','I_1','R_1','a_1','b_1','S_2','I_2','R_2','a_2','b_2','adj_1','adj_2','adj_3','adj_4','adj_5','adj_6','adj_7','adj_8','adj_9','S_3','I_3','R_3','a_3','b_3');