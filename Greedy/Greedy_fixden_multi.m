%Vaccine Allocation with Greedy Algorithm (Multiple network and fixed density) 
%Author: Guanyi Wang and Toru Kitagawa

clear all

% DATA INPUT AND PREPARATION
load('../simulation_multi_fixden.mat')

times = length(adj_1);%get the number of random networks
N_1 = length(S_1); 
N_2 = length(S_2);
N_3 = length(S_3);

%choose weight parameter
weight_1 = 1;
weight_2 = 1;

% Parameter setting 
parameter = true; % switch: choose 1st parameter set into account or 2nd parameter set
budget_constraint = 2; %swtich: choose the capacity constraint level from 1, 2, 3

if (budget_constraint == 1)
if (parameter)
diary('greedy_multi_fixden_para1_multi.log');
filename_results = '../greedy_multi_fixden_para1_bgt1.mat';
beta_11 = 0.7; %infection rate between group 1 and group 1
beta_12 = 0.5; %infection rate from group 2 to group 1
beta_21 = 0.5; %infection rate from group 1 to group 2
beta_22 = 0.6; %infection rate between group 2 and group 2
gamma_1 = 0.1; %recover rate for group 1
gamma_2 = 0.05; %recover rate for group 2
d_1 = round(0.07*N_1); %capacity constraint for N_1
d_2 = round(0.07*N_2); %capacity constraint for N_2
d_3 = round(0.07*N_3); %capacity constraint for N_3
else
diary('greedy_multi_fixden_para2.log');
filename_results = '../greedy_multi_fixden_para2_bgt1.mat';
beta_11 = 0.8; %infection rate between group 1 and group 1
beta_12 = 0.5; %infection rate from group 2 to group 1
beta_21 = 0.7; %infection rate from group 1 to group 2
beta_22 = 0.7; %infection rate between group 2 and group 2
gamma_1 = 0.1; %recover rate for group 1
gamma_2 = 0.025; %recover rate for group 2 
d_1 = round(0.07*N_1); %capacity constraint for N_1
d_2 = round(0.07*N_2); %capacity constraint for N_2
d_3 = round(0.07*N_3); %capacity constraint for N_3
end
end
if (budget_constraint == 2)
    if (parameter)
diary('greedy_multi_fixden_para1.log');
filename_results = '../greedy_multi_fixden_para1_bgt2.mat';
beta_11 = 0.7; %infection rate between group 1 and group 1
beta_12 = 0.5; %infection rate from group 2 to group 1
beta_21 = 0.5; %infection rate from group 1 to group 2
beta_22 = 0.6; %infection rate between group 2 and group 2
gamma_1 = 0.1; %recover rate for group 1
gamma_2 = 0.05; %recover rate for group 2
d_1 = round(0.1*N_1); %capacity constraint for N_1
d_2 = round(0.1*N_2); %capacity constraint for N_2
d_3 = round(0.1*N_3); %capacity constraint for N_3
else
diary('greedy_multi_fixden_para2.log');
filename_results = '../greedy_multi_fixden_para2_bgt2.mat';
beta_11 = 0.8; %infection rate between group 1 and group 1
beta_12 = 0.5; %infection rate from group 2 to group 1
beta_21 = 0.7; %infection rate from group 1 to group 2
beta_22 = 0.7; %infection rate between group 2 and group 2
gamma_1 = 0.1; %recover rate for group 1
gamma_2 = 0.025; %recover rate for group 2   
d_1 = round(0.1*N_1); %capacity constraint for N_1
d_2 = round(0.1*N_2); %capacity constraint for N_2
d_3 = round(0.1*N_3); %capacity constraint for N_3
    end
end

if (budget_constraint == 3)
    if (parameter)
diary('greedy_multi_fixden_para1.log');
filename_results = '../greedy_multi_fixden_para1_bgt3.mat';
beta_11 = 0.7; %infection rate between group 1 and group 1
beta_12 = 0.5; %infection rate from group 2 to group 1
beta_21 = 0.5; %infection rate from group 1 to group 2
beta_22 = 0.6; %infection rate between group 2 and group 2
gamma_1 = 0.1; %recover rate for group 1
gamma_2 = 0.05; %recover rate for group 2
d_1 = round(0.2*N_1); %capacity constraint for N_1
d_2 = round(0.2*N_2); %capacity constraint for N_2
d_3 = round(0.2*N_3); %capacity constraint for N_3
else
diary('greedy_multi_fixden_para2.log');
filename_results = '../greedy_multi_fixden_para2_bgt3.mat';
beta_11 = 0.8; %infection rate between group 1 and group 1
beta_12 = 0.5; %infection rate from group 2 to group 1
beta_21 = 0.7; %infection rate from group 1 to group 2
beta_22 = 0.7; %infection rate between group 2 and group 2
gamma_1 = 0.1; %recover rate for group 1
gamma_2 = 0.025; %recover rate for group 2  
d_1 = round(0.2*N_1); %capacity constraint for N_1
d_2 = round(0.2*N_2); %capacity constraint for N_2
d_3 = round(0.2*N_3); %capacity constraint for N_3
    end
end

%define the weight vector
Indx_1a  = (a_1 == 1);
Indx_1b  = (b_1 == 1);
g_1 = ones(N_1,1);%weight for N_1
g_1(Indx_1a) = g_1(Indx_1a)*weight_1;
g_1(Indx_1b) = g_1(Indx_1b)*weight_2;

Indx_2a  = (a_2 == 1);
Indx_2b  = (b_2 == 1);
g_2 = ones(N_2,1);%weight for N_2
g_2(Indx_2a) = g_2(Indx_2a)*weight_1;
g_2(Indx_2b) = g_2(Indx_2b)*weight_2;

Indx_3a = (a_3 == 1);
Indx_3b = (b_3 == 1);
g_3 = ones(N_3,1);%weight for N_3
g_3(Indx_3a) = g_3(Indx_3a)*weight_1;
g_3(Indx_3b) = g_3(Indx_3b)*weight_2;

tic %calculate the running time
for qq = 1:times
%Define the parameter vector
C_1 = g_1 - diag(R_1*g_1') - diag(diag(a_1*I_1')*g_1')*gamma_1 - diag(diag(b_1*I_1')*g_1')*gamma_2 - diag(S_1*g_1'); %define the vector for N_1
C_2 = g_2 - diag(R_2*g_2') - diag(diag(a_2*I_2')*g_2')*gamma_1 - diag(diag(b_2*I_2')*g_2')*gamma_2 - diag(S_2*g_2'); %define the vector for N_2
C_3 = g_3 - diag(R_3*g_3') - diag(diag(a_3*I_3')*g_3')*gamma_1 - diag(diag(b_3*I_3')*g_3')*gamma_2 - diag(S_3*g_3'); %define the vector for N_3
%Define the weighting matrix

%the first element vector in weighting matrix
%the first element vector in weighting matrix
W_11{qq} = beta_11*diag(g_1)*diag(a_1)*diag(S_1)*adj_1{qq}*diag(a_1)*diag(I_1);
W_12{qq} = beta_11*diag(g_1)*diag(a_1)*diag(S_1)*adj_2{qq}*diag(a_1)*diag(I_1);
W_13 = beta_11*diag(g_1)*diag(a_1)*diag(S_1)*adj_3*diag(a_1)*diag(I_1);
W_14{qq} = beta_11*diag(g_2)*diag(a_2)*diag(S_2)*adj_4{qq}*diag(a_2)*diag(I_2);
W_15{qq} = beta_11*diag(g_2)*diag(a_2)*diag(S_2)*adj_5{qq}*diag(a_2)*diag(I_2);
W_16 = beta_11*diag(g_2)*diag(a_2)*diag(S_2)*adj_6*diag(a_2)*diag(I_2);
W_17{qq} = beta_11*diag(g_3)*diag(a_3)*diag(S_3)*adj_7{qq}*diag(a_3)*diag(I_3);
W_18{qq} = beta_11*diag(g_3)*diag(a_3)*diag(S_3)*adj_8{qq}*diag(a_3)*diag(I_3);
W_19 = beta_11*diag(g_3)*diag(a_3)*diag(S_3)*adj_9*diag(a_3)*diag(I_3);

%the second element vector in weighting matrix
W_21{qq} = beta_12*diag(g_1)*diag(a_1)*diag(S_1)*adj_1{qq}*diag(b_1)*diag(I_1);
W_22{qq} = beta_12*diag(g_1)*diag(a_1)*diag(S_1)*adj_2{qq}*diag(b_1)*diag(I_1);
W_23 = beta_12*diag(g_1)*diag(a_1)*diag(S_1)*adj_3*diag(b_1)*diag(I_1);
W_24{qq} = beta_12*diag(g_2)*diag(a_2)*diag(S_2)*adj_4{qq}*diag(b_2)*diag(I_2);
W_25{qq} = beta_12*diag(g_2)*diag(a_2)*diag(S_2)*adj_5{qq}*diag(b_2)*diag(I_2);
W_26 = beta_12*diag(g_2)*diag(a_2)*diag(S_2)*adj_6*diag(b_2)*diag(I_2);
W_27{qq} = beta_12*diag(g_3)*diag(a_3)*diag(S_3)*adj_7{qq}*diag(b_3)*diag(I_3);
W_28{qq} = beta_12*diag(g_3)*diag(a_3)*diag(S_3)*adj_8{qq}*diag(b_3)*diag(I_3);
W_29 = beta_12*diag(g_3)*diag(a_3)*diag(S_3)*adj_9*diag(b_3)*diag(I_3);

%the third element vector in weighting matrix
W_31{qq} = beta_21*diag(g_1)*diag(b_1)*diag(S_1)*adj_1{qq}*diag(a_1)*diag(I_1);
W_32{qq} = beta_21*diag(g_1)*diag(b_1)*diag(S_1)*adj_2{qq}*diag(a_1)*diag(I_1);
W_33 = beta_21*diag(g_1)*diag(b_1)*diag(S_1)*adj_3*diag(a_1)*diag(I_1);
W_34{qq} = beta_21*diag(g_2)*diag(b_2)*diag(S_2)*adj_4{qq}*diag(a_2)*diag(I_2);
W_35{qq} = beta_21*diag(g_2)*diag(b_2)*diag(S_2)*adj_5{qq}*diag(a_2)*diag(I_2);
W_36 = beta_21*diag(g_2)*diag(b_2)*diag(S_2)*adj_6*diag(a_2)*diag(I_2);
W_37{qq} = beta_21*diag(g_3)*diag(b_3)*diag(S_3)*adj_7{qq}*diag(a_3)*diag(I_3);
W_38{qq} = beta_21*diag(g_3)*diag(b_3)*diag(S_3)*adj_8{qq}*diag(a_3)*diag(I_3);
W_39 = beta_21*diag(g_3)*diag(b_3)*diag(S_3)*adj_9*diag(a_3)*diag(I_3);

%the fourth element vector in weighting matrix
W_41{qq} = beta_22*diag(g_1)*diag(b_1)*diag(S_1)*adj_1{qq}*diag(b_1)*diag(I_1);
W_42{qq} = beta_22*diag(g_1)*diag(b_1)*diag(S_1)*adj_2{qq}*diag(b_1)*diag(I_1);
W_43 = beta_22*diag(g_1)*diag(b_1)*diag(S_1)*adj_3*diag(b_1)*diag(I_1);
W_44{qq} = beta_22*diag(g_2)*diag(b_2)*diag(S_2)*adj_4{qq}*diag(b_2)*diag(I_2);
W_45{qq} = beta_22*diag(g_2)*diag(b_2)*diag(S_2)*adj_5{qq}*diag(b_2)*diag(I_2);
W_46 = beta_22*diag(g_2)*diag(b_2)*diag(S_2)*adj_6*diag(b_2)*diag(I_2);
W_47{qq} = beta_22*diag(g_3)*diag(b_3)*diag(S_3)*adj_7{qq}*diag(b_3)*diag(I_3);
W_48{qq} = beta_22*diag(g_3)*diag(b_3)*diag(S_3)*adj_8{qq}*diag(b_3)*diag(I_3);
W_49 = beta_22*diag(g_3)*diag(b_3)*diag(S_3)*adj_9*diag(b_3)*diag(I_3);

%combine above elements together
W_1{qq} = - W_11{qq} - W_21{qq} - W_31{qq} - W_41{qq};
W_2{qq} = - W_12{qq} - W_22{qq} - W_32{qq} - W_42{qq};
W_3 = - W_13 - W_23 - W_33 - W_43;
W_4{qq} = - W_14{qq} - W_24{qq} - W_34{qq} - W_44{qq};
W_5{qq} = - W_15{qq} - W_25{qq} - W_35{qq} - W_45{qq};
W_6 = - W_16 - W_26 - W_36 - W_46;
W_7{qq} = - W_17{qq} - W_27{qq} - W_37{qq} - W_47{qq};
W_8{qq} = - W_18{qq} - W_28{qq} - W_38{qq} - W_48{qq};
W_9 = - W_19 - W_29 - W_39 - W_49;

%Normalization for the weighting matrix
for cc = 1:N_1
    if sum(adj_1{qq}(cc,:)) == 0 %for the unit with 0 neighbor, there is no normalization
       W_1{qq}(cc,:) = W_1{qq}(cc,:);
    else
       W_1{qq}(cc,:) = W_1{qq}(cc,:)/sum(adj_1{qq}(cc,:));
    end
    if sum(adj_2{qq}(cc,:)) == 0
        W_2{qq}(cc,:) = W_2{qq}(cc,:);
    else
        W_2{qq}(cc,:) = W_2{qq}(cc,:)/sum(adj_2{qq}(cc,:));
    end
 
    if sum(adj_3(cc,:)) == 0
        W_3(cc,:) = W_3(cc,:);
    else
        W_3(cc,:) = W_3(cc,:)/sum(adj_3(cc,:));
    end
end

for cc = 1:N_2
    if sum(adj_4{qq}(cc,:)) == 0
       W_4{qq}(cc,:) = W_4{qq}(cc,:);
    else
       W_4{qq}(cc,:) = W_4{qq}(cc,:)/sum(adj_4{qq}(cc,:));
    end
    if sum(adj_5{qq}(cc,:)) == 0
        W_5{qq}(cc,:) = W_5{qq}(cc,:);
    else
        W_5{qq}(cc,:) = W_5{qq}(cc,:)/sum(adj_5{qq}(cc,:));
    end
 
    if sum(adj_6(cc,:)) == 0
        W_6(cc,:) = W_6(cc,:);
    else
        W_6(cc,:) = W_6(cc,:)/sum(adj_6(cc,:));
    end
end

for cc = 1:N_3
    if sum(adj_7{qq}(cc,:)) == 0
       W_7{qq}(cc,:) = W_7{qq}(cc,:);
    else
       W_7{qq}(cc,:) = W_7{qq}(cc,:)/sum(adj_7{qq}(cc,:));
    end
    if sum(adj_8{qq}(cc,:)) == 0
        W_8{qq}(cc,:) = W_8{qq}(cc,:);
    else
        W_8{qq}(cc,:) = W_8{qq}(cc,:)/sum(adj_8{qq}(cc,:));
    end
 
    if sum(adj_9(cc,:)) == 0
        W_9(cc,:) = W_9(cc,:);
    else
        W_9(cc,:) = W_9(cc,:)/sum(adj_9(cc,:));
    end
end

%Define the index vector
Individual_11{qq} = (1:N_1)';
Individual_12{qq} = (1:N_1)';
Individual_13 = (1:N_1)';
Individual_21{qq} = (1:N_2)';
Individual_22{qq} = (1:N_2)';
Individual_23 = (1:N_2)';
Individual_31{qq} = (1:N_3)';
Individual_32{qq} = (1:N_3)';
Individual_33 = (1:N_3)';

%Define a vector to store the result
V_result_11{qq} = zeros(N_1,1);
V_result_12{qq} = zeros(N_1,1);
V_result_13 = zeros(N_1,1);
V_result_21{qq} = zeros(N_2,1);
V_result_22{qq} = zeros(N_2,1);
V_result_23 = zeros(N_2,1);
V_result_31{qq} = zeros(N_3,1);
V_result_32{qq} = zeros(N_3,1);
V_result_33 = zeros(N_3,1);

%Initialize a vector of V
V_11{qq} = zeros(N_1,d_1); 
V_12{qq} = zeros(N_1,d_1); 
V_13 = zeros(N_1,d_1); 
V_21{qq} = zeros(N_2,d_2); 
V_22{qq} = zeros(N_2,d_2); 
V_23 = zeros(N_2,d_2);
V_31{qq} = zeros(N_3,d_3); 
V_32{qq} = zeros(N_3,d_3); 
V_33 = zeros(N_3,d_3);

%Set a vector to store the outcome variable
F_11{qq} = zeros(N_1,d_1); 
F_12{qq} = zeros(N_1,d_1); 
F_13 = zeros(N_1,d_1); 
F_21{qq} = zeros(N_2,d_2); 
F_22{qq} = zeros(N_2,d_2); 
F_23 = zeros(N_2,d_2);
F_31{qq} = zeros(N_3,d_3); 
F_32{qq} = zeros(N_3,d_3); 
F_33 = zeros(N_3,d_3);

%set initial parameters for iteration
t_11 = 0;
g_11 = 1;
t_12 = 0;
g_12 = 1;
t_13 = 0;
g_13 = 1;
t_21 = 0;
g_21 = 1;
t_22 = 0;
g_22 = 1;
t_23 = 0;
g_23 = 1;
t_31 = 0;
g_31 = 1;
t_32 = 0;
g_32 = 1;
t_33 = 0;
g_33 = 1;

%start iteration for N_1, edges= 7N
while (g_11 <= d_1)
i = 1;
while (i < N_1-t_11+1)
    V_11{qq}(Individual_11{qq}(i),g_11) = 1;
    F_11{qq}(i,g_11) = obj(V_11{qq}(:,g_11),W_1{qq},C_1,g_1);
    V_11{qq}(Individual_11{qq}(i),g_11) = 0;
    i = i+1;
end

outcome_11{qq} = F_11{qq}(:,g_11);%generate a new vector of outcome in order to remove the redundent part

if (t_11>0)
    outcome_11{qq}(end-t_11+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_11{qq}] = sort(outcome_11{qq},'descend'); 
Individual_11{qq} = Individual_11{qq}(idx_11{qq});
V_result_11{qq}(g_11) = Individual_11{qq}(1); %record the vaccinated unit
Individual_11{qq}(1) = []; %drop the unit with highest marginal gain from ground set
Individual_11{qq} = sort(Individual_11{qq});
V_11{qq}(V_result_11{qq}(g_11),g_11:end) = 1;
t_11 = t_11+1;
g_11 = g_11+1;
end
%start iteration for N_1, edges= 10N
while (g_12 < d_1+1)
i = 1;
while (i < N_1-t_12+1)
    V_12{qq}(Individual_12{qq}(i),g_12) = 1;
    F_12{qq}(i,g_12) = obj(V_12{qq}(:,g_12),W_2{qq},C_1,g_1);
    V_12{qq}(Individual_12{qq}(i),g_12) = 0;
    i = i+1;
end

outcome_12{qq} = F_12{qq}(:,g_12);%generate a new vector of outcome in order to remove the redundent part

if (t_12>0)
    outcome_12{qq}(end-t_12+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_12{qq}] = sort(outcome_12{qq},'descend'); 
Individual_12{qq} = Individual_12{qq}(idx_12{qq});
V_result_12{qq}(g_12) = Individual_12{qq}(1); %record the vaccinated unit
Individual_12{qq}(1) = []; %drop the unit with highest marginal gain from ground set
Individual_12{qq} = sort(Individual_12{qq});
V_12{qq}(V_result_12{qq}(g_12),g_12:end) = 1;
t_12 = t_12+1;
g_12 = g_12+1;
end

%start iteration for N_1, edges= full
while (g_13 < d_1+1)
i = 1;
while (i < N_1-t_13+1)
    V_13(Individual_13(i),g_13) = 1;
    F_13(i,g_13) = obj(V_13(:,g_13),W_3,C_1,g_1);
    V_13(Individual_13(i),g_13) = 0;
    i = i+1;
end

outcome_13 = F_13(:,g_13);%generate a new vector of outcome in order to remove the redundent part

if (t_13>0)
    outcome_13(end-t_13+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_13] = sort(outcome_13,'descend'); 
Individual_13 = Individual_13(idx_13);
V_result_13(g_13) = Individual_13(1); %record the vaccinated unit
Individual_13(1) = []; %drop the unit with highest marginal gain from ground set
Individual_13 = sort(Individual_13);
V_13(V_result_13(g_13),g_13:end) = 1;
t_13 = t_13+1;
g_13 = g_13+1;
end

%start iteration for for N_2, edges= 7N
while (g_21 < d_2+1)
i = 1;
while (i < N_2-t_21+1)
    V_21{qq}(Individual_21{qq}(i),g_21) = 1;
    F_21{qq}(i,g_21) = obj(V_21{qq}(:,g_21),W_4{qq},C_2,g_2);
    V_21{qq}(Individual_21{qq}(i),g_21) = 0;
    i = i+1;
end

outcome_21{qq} = F_21{qq}(:,g_21);%generate a new vector of outcome in order to remove the redundent part

if (t_21>0)
    outcome_21{qq}(end-t_21+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_21{qq}] = sort(outcome_21{qq},'descend'); 
Individual_21{qq} = Individual_21{qq}(idx_21{qq});
V_result_21{qq}(g_21) = Individual_21{qq}(1); %record the vaccinated unit
Individual_21{qq}(1) = []; %drop the unit with highest marginal gain from ground set
Individual_21{qq} = sort(Individual_21{qq});
V_21{qq}(V_result_21{qq}(g_21),g_21:end) = 1;
t_21 = t_21+1;
g_21 = g_21+1;
end

%start iteration for N_2, edges= 10N
while (g_22 < d_2+1)
i = 1;
while (i < N_2-t_22+1)
    V_22{qq}(Individual_22{qq}(i),g_22) = 1;
    F_22{qq}(i,g_22) = obj(V_22{qq}(:,g_22),W_5{qq},C_2,g_2);
    V_22{qq}(Individual_22{qq}(i),g_22) = 0;
    i = i+1;
end

outcome_22{qq} = F_22{qq}(:,g_22);%generate a new vector of outcome in order to remove the redundent part

if (t_22>0)
    outcome_22{qq}(end-t_22+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_22{qq}] = sort(outcome_22{qq},'descend'); 
Individual_22{qq} = Individual_22{qq}(idx_22{qq});
V_result_22{qq}(g_22) = Individual_22{qq}(1); %record the vaccinated unit
Individual_22{qq}(1) = []; %drop the unit with highest marginal gain from ground set
Individual_22{qq} = sort(Individual_22{qq});
V_22{qq}(V_result_22{qq}(g_22),g_22:end) = 1;
t_22 = t_22+1;
g_22 = g_22+1;
end
%start iteration for N_2, edges= full
while (g_23 < d_2+1)
i = 1;
while (i < N_2-t_23+1)
    V_23(Individual_23(i),g_23) = 1;
    F_23(i,g_23) = obj(V_23(:,g_23),W_6,C_2,g_2);
    V_23(Individual_23(i),g_23) = 0;
    i = i+1;
end

outcome_23 = F_23(:,g_23);%generate a new vector of outcome in order to remove the redundent part

if (t_23>0)
    outcome_23(end-t_23+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_23] = sort(outcome_23,'descend'); 
Individual_23 = Individual_23(idx_23);
V_result_23(g_23) = Individual_23(1); %record the vaccinated unit
Individual_23(1) = []; %drop the unit with highest marginal gain from ground set
Individual_23 = sort(Individual_23);
V_23(V_result_23(g_23),g_23:end) = 1;
t_23 = t_23+1;
g_23 = g_23+1;
end

%start iteration for N_3, edges= 7N
while (g_31 < d_3+1)
i = 1;
while (i < N_3-t_31+1)
    V_31{qq}(Individual_31{qq}(i),g_31) = 1;
    F_31{qq}(i,g_31) = obj(V_31{qq}(:,g_31),W_7{qq},C_3,g_3);
    V_31{qq}(Individual_31{qq}(i),g_31) = 0;
    i = i+1;
end

outcome_31{qq} = F_31{qq}(:,g_31);%generate a new vector of outcome in order to remove the redundent part

if (t_31>0)
    outcome_31{qq}(end-t_31+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_31{qq}] = sort(outcome_31{qq},'descend'); 
Individual_31{qq} = Individual_31{qq}(idx_31{qq});
V_result_31{qq}(g_31) = Individual_31{qq}(1); %record the vaccinated unit
Individual_31{qq}(1) = []; %drop the unit with highest marginal gain from ground set
Individual_31{qq} = sort(Individual_31{qq});
V_31{qq}(V_result_31{qq}(g_31),g_31:end) = 1;
t_31 = t_31+1;
g_31 = g_31+1;
end

%start iteration for N_3, edges= 10N
while (g_32 < d_3+1)
i = 1;
while (i < N_3-t_32+1)
    V_32{qq}(Individual_32{qq}(i),g_32) = 1;
    F_32{qq}(i,g_32) = obj(V_32{qq}(:,g_32),W_8{qq},C_3,g_3);
    V_32{qq}(Individual_32{qq}(i),g_32) = 0;
    i = i+1;
end

outcome_32{qq} = F_32{qq}(:,g_32);%generate a new vector of outcome in order to remove the redundent part

if (t_32>0)
    outcome_32{qq}(end-t_32+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_32{qq}] = sort(outcome_32{qq},'descend'); 
Individual_32{qq} = Individual_32{qq}(idx_32{qq});
V_result_32{qq}(g_32) = Individual_32{qq}(1); %record the vaccinated unit
Individual_32{qq}(1) = []; %drop the unit with highest marginal gain from ground set
Individual_32{qq} = sort(Individual_32{qq});
V_32{qq}(V_result_32{qq}(g_32),g_32:end) = 1;
t_32 = t_32+1;
g_32 = g_32+1;
end

%start iteration for N_3, edges= full
while (g_33 < d_3+1)
i = 1;
while (i < N_3-t_33+1)
    V_33(Individual_33(i),g_33) = 1;
    F_33(i,g_33) = obj(V_33(:,g_33),W_9,C_3,g_3);
    V_33(Individual_33(i),g_33) = 0;
    i = i+1;
end

outcome_33 = F_33(:,g_33);%generate a new vector of outcome in order to remove the redundent part

if (t_33>0)
    outcome_33(end-t_33+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_33] = sort(outcome_33,'descend'); 
Individual_33 = Individual_33(idx_33);
V_result_33(g_33) = Individual_33(1); %record the vaccinated unit
Individual_33(1) = []; %drop the unit with highest marginal gain from ground set
Individual_33 = sort(Individual_33);
V_33(V_result_33(g_33),g_33:end) = 1;
t_33 = t_33+1;
g_33 = g_33+1;
end

%remove zero element from result
V_result_11{qq} = nonzeros(V_result_11{qq});
V_result_12{qq} = nonzeros(V_result_12{qq});
V_result_13 = nonzeros(V_result_13);
V_result_21{qq} = nonzeros(V_result_21{qq});
V_result_22{qq} = nonzeros(V_result_22{qq});
V_result_23 = nonzeros(V_result_23);
V_result_31{qq} = nonzeros(V_result_31{qq});
V_result_32{qq} = nonzeros(V_result_32{qq});
V_result_33 = nonzeros(V_result_33);
%sort the vaccinated unit
V_result_11{qq} = sort(V_result_11{qq});
V_result_12{qq} = sort(V_result_12{qq});
V_result_13 = sort(V_result_13);
V_result_21{qq} = sort(V_result_21{qq});
V_result_22{qq} = sort(V_result_22{qq});
V_result_23 = sort(V_result_23);
V_result_31{qq} = sort(V_result_31{qq});
V_result_32{qq} = sort(V_result_32{qq});
V_result_33 = sort(V_result_33);
%save the final result of social security
F_result_11{qq} = max(F_11{qq}(:,d_1));
F_result_12{qq} = max(F_12{qq}(:,d_1));
F_result_13 = max(F_13(:,d_1));
F_result_21{qq} = max(F_21{qq}(:,d_2));
F_result_22{qq} = max(F_22{qq}(:,d_2));
F_result_23 = max(F_23(:,d_2));
F_result_31{qq} = max(F_31{qq}(:,d_3));
F_result_32{qq} = max(F_32{qq}(:,d_3));
F_result_33 = max(F_33(:,d_3));

%calculate the constant part in the initial objective function
%constant part without interaction, N_1
zz = 1;
const_1 = 0;
while (zz<=N_1)
   const_1 = const_1 + R_1(zz)+a_1(zz)*I_1(zz)*gamma_1+b_1(zz)*I_1(zz)*gamma_2+S_1(zz);
    zz = zz+1;
end
%take average for the constant part without interaction 
const_1 = const_1/N_1;

%calcuate the interacted constant part sepearately, N_1 
mm = 1;
nn = 1;
const_11{qq} = 0;
while (nn <= N_1)
while (mm <= N_1)
   const_11{qq} = const_11{qq}+(beta_11*a_1(nn)*a_1(mm)+beta_12*a_1(nn)*b_1(mm)+beta_21*b_1(nn)*a_1(mm)+ beta_22*b_1(nn)*b_1(mm))*adj_1{qq}(nn,mm)*I_1(mm)*S_1(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_1
const_11{qq} = const_11{qq}/N_1;

mm = 1;
nn = 1;
const_12{qq} = 0;
while (nn <= N_1)
while (mm <= N_1)
   const_12{qq} = const_12{qq}+(beta_11*a_1(nn)*a_1(mm)+beta_12*a_1(nn)*b_1(mm)+beta_21*b_1(nn)*a_1(mm)+ beta_22*b_1(nn)*b_1(mm))*adj_2{qq}(nn,mm)*I_1(mm)*S_1(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_2
const_12{qq} = const_12{qq}/N_1;

mm = 1;
nn = 1;
const_13 = 0;
while (nn <= N_1)
while (mm <= N_1)
   const_13 = const_13+(beta_11*a_1(nn)*a_1(mm)+beta_12*a_1(nn)*b_1(mm)+beta_21*b_1(nn)*a_1(mm)+ beta_22*b_1(nn)*b_1(mm))*adj_3(nn,mm)*I_1(mm)*S_1(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_3
const_13 = const_13/N_1;

%constant part without interaction, N_2
zz = 1;
const_2 = 0;
while (zz<=N_2)
   const_2 = const_2 + R_2(zz)+a_2(zz)*I_2(zz)*gamma_2+b_2(zz)*I_2(zz)*gamma_2+S_2(zz);
    zz = zz+1;
end
%take average for the constant part without interaction 
const_2 = const_2/N_2;

%calcuate the interacted constant part sepearately, N_2
mm = 1;
nn = 1;
const_21{qq} = 0;
while (nn <= N_2)
while (mm <= N_2)
   const_21{qq} = const_21{qq}+(beta_11*a_2(nn)*a_2(mm)+beta_12*a_2(nn)*b_2(mm)+beta_21*b_2(nn)*a_2(mm)+ beta_22*b_2(nn)*b_2(mm))*adj_4{qq}(nn,mm)*I_2(mm)*S_2(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_4
const_21{qq} = const_21{qq}/N_2;

mm = 1;
nn = 1;
const_22{qq} = 0;
while (nn <= N_2)
while (mm <= N_2)
   const_22{qq} = const_22{qq} +(beta_11*a_2(nn)*a_2(mm)+beta_12*a_2(nn)*b_2(mm)+beta_21*b_2(nn)*a_2(mm)+ beta_22*b_2(nn)*b_2(mm))*adj_5{qq}(nn,mm)*I_2(mm)*S_2(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_5
const_22{qq} = const_22{qq}/N_2;

mm = 1;
nn = 1;
const_23 = 0;
while (nn <= N_2)
while (mm <= N_2)
   const_23 = const_23+(beta_11*a_2(nn)*a_2(mm)+beta_12*a_2(nn)*b_2(mm)+beta_21*b_2(nn)*a_2(mm)+ beta_22*b_2(nn)*b_2(mm))*adj_6(nn,mm)*I_2(mm)*S_2(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_6
const_23 = const_23/N_2;

%constant part without interaction, N_3
zz = 1;
const_3 = 0;
while (zz<=N_3)
   const_3 = const_3 + R_3(zz)+a_3(zz)*I_3(zz)*gamma_1+b_3(zz)*I_3(zz)*gamma_2+S_3(zz);
    zz = zz+1;
end

%take average for the constant part without interaction 
const_3 = const_3/N_3;

%calcuate the interacted constant part sepearately, N_3
mm = 1;
nn = 1;
const_31{qq} = 0;
while (nn <= N_3)
while (mm <= N_3)
   const_31{qq} = const_31{qq}+(beta_11*a_3(nn)*a_3(mm)+beta_12*a_3(nn)*b_3(mm)+beta_21*b_3(nn)*a_3(mm)+ beta_22*b_3(nn)*b_3(mm))*adj_7{qq}(nn,mm)*I_3(mm)*S_3(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_7
const_31{qq} = const_31{qq}/N_3;

mm = 1;
nn = 1;
const_32{qq} = 0;
while (nn <= N_3)
while (mm <= N_3)
   const_32{qq} = const_32{qq}+(beta_11*a_3(nn)*a_3(mm)+beta_12*a_3(nn)*b_3(mm)+beta_21*b_3(nn)*a_3(mm)+ beta_22*b_3(nn)*b_3(mm))*adj_8{qq}(nn,mm)*I_3(mm)*S_3(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_8
const_32{qq} = const_32{qq}/N_3;

mm = 1;
nn = 1;
const_33 = 0;
while (nn <= N_3)
while (mm <= N_3)
   const_33 = const_33+(beta_11*a_3(nn)*a_3(mm)+beta_12*a_3(nn)*b_3(mm)+beta_21*b_3(nn)*a_3(mm)+ beta_22*b_3(nn)*b_3(mm))*adj_9(nn,mm)*I_3(mm)*S_3(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_9
const_33 = const_33/N_3;
%combime the previous result together
Y_11{qq} = const_1-const_11{qq}+F_result_11{qq};
Y_12{qq} = const_1-const_12{qq}+F_result_12{qq};
Y_13 = const_1-const_13+F_result_13;
Y_21{qq} = const_2-const_21{qq}+F_result_21{qq};
Y_22{qq} = const_2-const_22{qq}+F_result_22{qq};
Y_23 = const_2-const_23+F_result_23;
Y_31{qq} = const_3-const_31{qq}+F_result_31{qq};
Y_32{qq} = const_3-const_32{qq}+F_result_32{qq};
Y_33 = const_3-const_33+F_result_33;
end

%calculate the average value of outcome variable for multiple network
%set initial value for the interation 
Y_11_ave = 0;
Y_12_ave = 0;
Y_13_ave = Y_13;
Y_21_ave = 0;
Y_22_ave = 0;
Y_23_ave = Y_23;
Y_31_ave = 0;
Y_32_ave = 0;
Y_33_ave = Y_33;

%take average value across different networks
ww = 1;
while (ww<=times)
Y_11_ave = Y_11_ave + Y_11{ww};
ww = ww+1;
end
Y_11_ave = Y_11_ave/times;

ww = 1;
while (ww<=times)
Y_12_ave = Y_12_ave + Y_12{ww};
ww = ww+1;
end
Y_12_ave = Y_12_ave/times;

ww = 1;
while (ww<=times)
Y_21_ave = Y_21_ave + Y_21{ww};
ww = ww+1;
end
Y_21_ave = Y_21_ave/times;

ww = 1;
while (ww<=times)
Y_22_ave = Y_22_ave + Y_22{ww};
ww = ww+1;
end
Y_22_ave = Y_22_ave/times;

ww = 1;
while (ww<=times)
Y_31_ave = Y_31_ave + Y_31{ww};
ww = ww+1;
end
Y_31_ave = Y_31_ave/times;

ww = 1;
while (ww<=times)
Y_32_ave = Y_32_ave + Y_32{ww};
ww = ww+1;
end
Y_32_ave = Y_32_ave/times;

%calculate standard error from greedy alogrithm on multiple networks
Y_11_dif = zeros(times,1);
for qq = 1:times
   Y_11_dif(qq) = (Y_11{qq} - Y_11_ave)^2;
end
Y_11_std = sqrt(sum(Y_11_dif));%standard error

Y_12_dif = zeros(times,1); 
for qq = 1:times
   Y_12_dif(qq) = (Y_12{qq} - Y_12_ave)^2;
end
Y_12_std = sqrt(sum(Y_12_dif));%standard error

Y_21_dif = zeros(times,1);
for qq = 1:times
   Y_21_dif(qq) = (Y_21{qq} - Y_21_ave)^2;
end
Y_21_std = sqrt(sum(Y_21_dif));%standard error

Y_22_dif = zeros(times,1);
for qq = 1:times
   Y_22_dif(qq) = (Y_22{qq} - Y_22_ave)^2;
end
Y_22_std = sqrt(sum(Y_22_dif));%standard error

Y_31_dif = zeros(times,1);
for qq = 1:times
   Y_31_dif(qq) = (Y_31{qq} - Y_31_ave)^2;
end
Y_31_std = sqrt(sum(Y_31_dif));%standard error

Y_32_dif = zeros(times,1);
for qq = 1:times
   Y_32_dif(qq) = (Y_32{qq} - Y_32_ave)^2;
end
Y_32_std = sqrt(sum(Y_32_dif));%standard error
toc
%clear the irrelevant variables
clear Individual_11 Individual_12 Individual_13 Individual_21 Individual_22 Individual_23 Individual_31 Individual_32 Individual_33
clear outcome_11 outcome_12 outcome_13 outcome_21 outcome_22 outcome_23 outcome_31 outcome_32 outcome_33
clear idx_11 idx_12 idx_13 idx_21 idx_22 idx_23 idx_31 idx_32 idx_33
clear t_11 t_12 t_13 t_21 t_22 t_23 t_31 t_32 t_33
clear i mm nn zz ww qq 
clear g_11 g_12 g_13 g_21 g_22 g_23 g_31 g_32 g_33
clear V_11 V_12 V_13 V_21 V_22 V_23 V_31 V_32 V_33
clear W_11 W_12 W_13 W_14 W_15 W_16 W_17 W_18 W_19 W_21 W_22 W_23 W_24 W_25 W_26 W_27 W_28 W_29 
clear W_31 W_32 W_33 W_34 W_35 W_36 W_37 W_38 W_39 W_41 W_42 W_43 W_44 W_45 W_46 W_47 W_48 W_49 
clear W_1 W_2 W_3 W_4 W_5 W_6 W_7 W_8 W_9
clear beta_11 beta_12 beta_21 beta_22 gamma_1 gamma_2
clear a_1 a_2 a_3 b_1 b_2 b_3 S_1 S_2 S_3 I_1 I_2 I_3 R_1 R_2 R_3
clear const_1 const_2 const_3 const_11 const_12 const_13 const_21 const_22 const_23 const_31 const_32 const_33
clear F_result_11 F_result_12 F_result_13 F_result_21 F_result_22 F_result_23 F_result_31 F_result_32 F_result_33
clear F_11 F_12 F_13 F_21 F_22 F_23 F_31 F_32 F_33 C_1 C_2 C_3 Y_11 Y_12 Y_13 Y_21 Y_22 Y_23 Y_31 Y_32 Y_33
clear G_1 G_2 G_3 G_4 G_5 G_6 G_7 G_8 G_9 adj_1 adj_2 adj_3 adj_4 adj_5 adj_6 adj_7 adj_8 adj_9
clear Y_11_dif Y_12_dif Y_21_dif Y_22_dif Y_31_dif Y_32_dif
% Save coefficients and other variables for graphs
save(filename_results);

diary off