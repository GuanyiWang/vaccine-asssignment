%Vaccine Allocation with Greedy Algorithm
%Author: Guanyi Wang and Toru Kitagawa

clear all


% DATA INPUT AND PREPARATION
load('../simulation.mat')

N_1 = length(S_1); 
N_2 = length(S_2);
N_3 = length(S_3);

% Parameter setting 
parameter = false; % switch: choose 1st parameter set into account or 2nd parameter set
budget_constraint = 2; %swtich: choose the capacity constraint level from 1, 2, 3

if (budget_constraint == 1)
if (parameter)
diary('greedy_para1.log');
filename_results = '../greedy_para1_bgt1.mat';
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
diary('greedy_para2.log');
filename_results = '../greedy_para2_bgt1.mat';
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
diary('greedy_para1.log');
filename_results = '../greedy_para1_bgt2.mat';
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
diary('greedy_para2.log');
filename_results = '../greedy_para2_bgt2.mat';
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
diary('greedy_para1.log');
filename_results = '../greedy_para1_bgt3.mat';
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
diary('greedy_para2.log');
filename_results = '../greedy_para2_bgt3.mat';
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

%Define the parameter vector
C_1 = ones(N_1,1) - R_1 - diag(a_1*I_1')*gamma_1 - diag(b_1*I_1')*gamma_2 - S_1; %define the vector for N_1
C_2 = ones(N_2,1) - R_2 - diag(a_2*I_2')*gamma_1 - diag(b_2*I_2')*gamma_2 - S_2; %define the vector for N_2
C_3 = ones(N_3,1) - R_3 - diag(a_3*I_3')*gamma_1 - diag(b_3*I_3')*gamma_2 - S_3; %define the vector for N_3

%Define the weighting matrix

%the first element vector in weighting matrix
W_11 = beta_11*diag(a_1)*diag(S_1)*adj_1*diag(a_1)*diag(I_1);
W_12 = beta_11*diag(a_1)*diag(S_1)*adj_2*diag(a_1)*diag(I_1);
W_13 = beta_11*diag(a_1)*diag(S_1)*adj_3*diag(a_1)*diag(I_1);
W_14 = beta_11*diag(a_2)*diag(S_2)*adj_4*diag(a_2)*diag(I_2);
W_15 = beta_11*diag(a_2)*diag(S_2)*adj_5*diag(a_2)*diag(I_2);
W_16 = beta_11*diag(a_2)*diag(S_2)*adj_6*diag(a_2)*diag(I_2);
W_17 = beta_11*diag(a_3)*diag(S_3)*adj_7*diag(a_3)*diag(I_3);
W_18 = beta_11*diag(a_3)*diag(S_3)*adj_8*diag(a_3)*diag(I_3);
W_19 = beta_11*diag(a_3)*diag(S_3)*adj_9*diag(a_3)*diag(I_3);

%the second element vector in weighting matrix
W_21 = beta_12*diag(a_1)*diag(S_1)*adj_1*diag(b_1)*diag(I_1);
W_22 = beta_12*diag(a_1)*diag(S_1)*adj_2*diag(b_1)*diag(I_1);
W_23 = beta_12*diag(a_1)*diag(S_1)*adj_3*diag(b_1)*diag(I_1);
W_24 = beta_12*diag(a_2)*diag(S_2)*adj_4*diag(b_2)*diag(I_2);
W_25 = beta_12*diag(a_2)*diag(S_2)*adj_5*diag(b_2)*diag(I_2);
W_26 = beta_12*diag(a_2)*diag(S_2)*adj_6*diag(b_2)*diag(I_2);
W_27 = beta_12*diag(a_3)*diag(S_3)*adj_7*diag(b_3)*diag(I_3);
W_28 = beta_12*diag(a_3)*diag(S_3)*adj_8*diag(b_3)*diag(I_3);
W_29 = beta_12*diag(a_3)*diag(S_3)*adj_9*diag(b_3)*diag(I_3);

%the third element vector in weighting matrix
W_31 = beta_21*diag(b_1)*diag(S_1)*adj_1*diag(a_1)*diag(I_1);
W_32 = beta_21*diag(b_1)*diag(S_1)*adj_2*diag(a_1)*diag(I_1);
W_33 = beta_21*diag(b_1)*diag(S_1)*adj_3*diag(a_1)*diag(I_1);
W_34 = beta_21*diag(b_2)*diag(S_2)*adj_4*diag(a_2)*diag(I_2);
W_35 = beta_21*diag(b_2)*diag(S_2)*adj_5*diag(a_2)*diag(I_2);
W_36 = beta_21*diag(b_2)*diag(S_2)*adj_6*diag(a_2)*diag(I_2);
W_37 = beta_21*diag(b_3)*diag(S_3)*adj_7*diag(a_3)*diag(I_3);
W_38 = beta_21*diag(b_3)*diag(S_3)*adj_8*diag(a_3)*diag(I_3);
W_39 = beta_21*diag(b_3)*diag(S_3)*adj_9*diag(a_3)*diag(I_3);

%the fourth element vector in weighting matrix
W_41 = beta_22*diag(b_1)*diag(S_1)*adj_1*diag(b_1)*diag(I_1);
W_42 = beta_22*diag(b_1)*diag(S_1)*adj_2*diag(b_1)*diag(I_1);
W_43 = beta_22*diag(b_1)*diag(S_1)*adj_3*diag(b_1)*diag(I_1);
W_44 = beta_22*diag(b_2)*diag(S_2)*adj_4*diag(b_2)*diag(I_2);
W_45 = beta_22*diag(b_2)*diag(S_2)*adj_5*diag(b_2)*diag(I_2);
W_46 = beta_22*diag(b_2)*diag(S_2)*adj_6*diag(b_2)*diag(I_2);
W_47 = beta_22*diag(b_3)*diag(S_3)*adj_7*diag(b_3)*diag(I_3);
W_48 = beta_22*diag(b_3)*diag(S_3)*adj_8*diag(b_3)*diag(I_3);
W_49 = beta_22*diag(b_3)*diag(S_3)*adj_9*diag(b_3)*diag(I_3);

%combine above elements together
W_1 = - W_11 - W_21 - W_31 - W_41;
W_2 = - W_12 - W_22 - W_32 - W_42;
W_3 = - W_13 - W_23 - W_33 - W_43;
W_4 = - W_14 - W_24 - W_34 - W_44;
W_5 = - W_15 - W_25 - W_35 - W_45;
W_6 = - W_16 - W_26 - W_36 - W_46;
W_7 = - W_17 - W_27 - W_37 - W_47;
W_8 = - W_18 - W_28 - W_38 - W_48;
W_9 = - W_19 - W_29 - W_39 - W_49;

%Normalization for the weighting matrix
for cc = 1:N_1
    if sum(adj_1(cc,:)) == 0 %for the unit with 0 neighbor, there is no normalization
       W_1(cc,:) = W_1(cc,:);
    else
       W_1(cc,:) = W_1(cc,:)/sum(adj_1(cc,:));
    end
    if sum(adj_2(cc,:)) == 0
        W_2(cc,:) = W_2(cc,:);
    else
        W_2(cc,:) = W_2(cc,:)/sum(adj_2(cc,:));
    end
 
    if sum(adj_3(cc,:)) == 0
        W_3(cc,:) = W_3(cc,:);
    else
        W_3(cc,:) = W_3(cc,:)/sum(adj_3(cc,:));
    end
end

for cc = 1:N_2
    if sum(adj_4(cc,:)) == 0
       W_4(cc,:) = W_4(cc,:);
    else
       W_4(cc,:) = W_4(cc,:)/sum(adj_4(cc,:));
    end
    if sum(adj_5(cc,:)) == 0
        W_5(cc,:) = W_5(cc,:);
    else
        W_5(cc,:) = W_5(cc,:)/sum(adj_5(cc,:));
    end
 
    if sum(adj_6(cc,:)) == 0
        W_6(cc,:) = W_6(cc,:);
    else
        W_6(cc,:) = W_6(cc,:)/sum(adj_6(cc,:));
    end
end

for cc = 1:N_3
    if sum(adj_7(cc,:)) == 0
       W_7(cc,:) = W_7(cc,:);
    else
       W_7(cc,:) = W_7(cc,:)/sum(adj_7(cc,:));
    end
    if sum(adj_8(cc,:)) == 0
        W_8(cc,:) = W_8(cc,:);
    else
        W_8(cc,:) = W_8(cc,:)/sum(adj_8(cc,:));
    end
 
    if sum(adj_9(cc,:)) == 0
        W_9(cc,:) = W_9(cc,:);
    else
        W_9(cc,:) = W_9(cc,:)/sum(adj_9(cc,:));
    end
end

%define the index vector
Individual_11 = (1:N_1)';
Individual_12 = (1:N_1)';
Individual_13 = (1:N_1)';
Individual_21 = (1:N_2)';
Individual_22 = (1:N_2)';
Individual_23 = (1:N_2)';
Individual_31 = (1:N_3)';
Individual_32 = (1:N_3)';
Individual_33 = (1:N_3)';
%vector for store the result
V_result_11 = zeros(N_1,1);
V_result_12 = zeros(N_1,1);
V_result_13 = zeros(N_1,1);
V_result_21 = zeros(N_2,1);
V_result_22 = zeros(N_2,1);
V_result_23 = zeros(N_2,1);
V_result_31 = zeros(N_3,1);
V_result_32 = zeros(N_3,1);
V_result_33 = zeros(N_3,1);
%Initialize the vector of V
V_11 = zeros(N_1,d_1); 
V_12 = zeros(N_1,d_1); 
V_13 = zeros(N_1,d_1); 
V_21 = zeros(N_2,d_2); 
V_22 = zeros(N_2,d_2); 
V_23 = zeros(N_2,d_2);
V_31 = zeros(N_3,d_3); 
V_32 = zeros(N_3,d_3); 
V_33 = zeros(N_3,d_3);
%Set the vector of objective function
F_11 = zeros(N_1,d_1); 
F_12 = zeros(N_1,d_1); 
F_13 = zeros(N_1,d_1); 
F_21 = zeros(N_2,d_2); 
F_22 = zeros(N_2,d_2); 
F_23 = zeros(N_2,d_2);
F_31 = zeros(N_3,d_3); 
F_32 = zeros(N_3,d_3); 
F_33 = zeros(N_3,d_3);
tic %calculate the running time
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
%start iteration for N_1, edges = 7N
while (g_11 <= d_1)
i = 1;
while (i < N_1-t_11+1)
    V_11(Individual_11(i),g_11) = 1;
    F_11(i,g_11) = obj(V_11(:,g_11),W_1,C_1);
    V_11(Individual_11(i),g_11) = 0;
    i = i+1;
end

outcome_11 = F_11(:,g_11);%generate a new vector of outcome in order to remove the redundent part

if (t_11>0)
    outcome_11(end-t_11+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_11] = sort(outcome_11,'descend'); 
Individual_11 = Individual_11(idx_11);
V_result_11(g_11) = Individual_11(1); %record the vaccinated unit
Individual_11(1) = []; %drop the unit with highest marginal gain from ground set
Individual_11 = sort(Individual_11);
V_11(V_result_11(g_11),g_11:end) = 1;
t_11 = t_11+1;
g_11 = g_11+1;
end
%start iteration for N_1, edges = 10N
while (g_12 < d_1+1)
i = 1;
while (i < N_1-t_12+1)
    V_12(Individual_12(i),g_12) = 1;
    F_12(i,g_12) = obj(V_12(:,g_12),W_2,C_1);
    V_12(Individual_12(i),g_12) = 0;
    i = i+1;
end

outcome_12 = F_12(:,g_12);%generate a new vector of outcome in order to remove the redundent part

if (t_12>0)
    outcome_12(end-t_12+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_12] = sort(outcome_12,'descend'); 
Individual_12 = Individual_12(idx_12);
V_result_12(g_12) = Individual_12(1); %record the vaccinated unit
Individual_12(1) = []; %drop the unit with highest marginal gain from ground set
Individual_12 = sort(Individual_12);
V_12(V_result_12(g_12),g_12:end) = 1;
t_12 = t_12+1;
g_12 = g_12+1;
end
%start iteration for N_1, edges = full
while (g_13 < d_1+1)
i = 1;
while (i < N_1-t_13+1)
    V_13(Individual_13(i),g_13) = 1;
    F_13(i,g_13) = obj(V_13(:,g_13),W_3,C_1);
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

%start iteration for N_2, edges = 7N
while (g_21 < d_2+1)
i = 1;
while (i < N_2-t_21+1)
    V_21(Individual_21(i),g_21) = 1;
    F_21(i,g_21) = obj(V_21(:,g_21),W_4,C_2);
    V_21(Individual_21(i),g_21) = 0;
    i = i+1;
end

outcome_21 = F_21(:,g_21);%generate a new vector of outcome in order to remove the redundent part

if (t_21>0)
    outcome_21(end-t_21+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_21] = sort(outcome_21,'descend'); 
Individual_21 = Individual_21(idx_21);
V_result_21(g_21) = Individual_21(1); %record the vaccinated unit
Individual_21(1) = []; %drop the unit with highest marginal gain from ground set
Individual_21 = sort(Individual_21);
V_21(V_result_21(g_21),g_21:end) = 1;
t_21 = t_21+1;
g_21 = g_21+1;
end
%start iteration for N_2, edges = 10N
while (g_22 < d_2+1)
i = 1;
while (i < N_2-t_22+1)
    V_22(Individual_22(i),g_22) = 1;
    F_22(i,g_22) = obj(V_22(:,g_22),W_5,C_2);
    V_22(Individual_22(i),g_22) = 0;
    i = i+1;
end

outcome_22 = F_22(:,g_22);%generate a new vector of outcome in order to remove the redundent part

if (t_22>0)
    outcome_22(end-t_22+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_22] = sort(outcome_22,'descend'); 
Individual_22 = Individual_22(idx_22);
V_result_22(g_22) = Individual_22(1); %record the vaccinated unit
Individual_22(1) = []; %drop the unit with highest marginal gain from ground set
Individual_22 = sort(Individual_22);
V_22(V_result_22(g_22),g_22:end) = 1;
t_22 = t_22+1;
g_22 = g_22+1;
end
%start iteration for N_2, edges = full
while (g_23 < d_2+1)
i = 1;
while (i < N_2-t_23+1)
    V_23(Individual_23(i),g_23) = 1;
    F_23(i,g_23) = obj(V_23(:,g_23),W_6,C_2);
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

%start iteration for N_3, edges = 7N
while (g_31 < d_3+1)
i = 1;
while (i < N_3-t_31+1)
    V_31(Individual_31(i),g_31) = 1;
    F_31(i,g_31) = obj(V_31(:,g_31),W_7,C_3);
    V_31(Individual_31(i),g_31) = 0;
    i = i+1;
end

outcome_31 = F_31(:,g_31);%generate a new vector of outcome in order to remove the redundent part

if (t_31>0)
    outcome_31(end-t_31+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_31] = sort(outcome_31,'descend'); 
Individual_31 = Individual_31(idx_31);
V_result_31(g_31) = Individual_31(1); %record the vaccinated unit
Individual_31(1) = []; %drop the unit with highest marginal gain from ground set
Individual_31 = sort(Individual_31);
V_31(V_result_31(g_31),g_31:end) = 1;
t_31 = t_31+1;
g_31 = g_31+1;
end

%start iteration for N_3, edges = 10N
while (g_32 < d_3+1)
i = 1;
while (i < N_3-t_32+1)
    V_32(Individual_32(i),g_32) = 1;
    F_32(i,g_32) = obj(V_32(:,g_32),W_8,C_3);
    V_32(Individual_32(i),g_32) = 0;
    i = i+1;
end

outcome_32 = F_32(:,g_32);%generate a new vector of outcome in order to remove the redundent part

if (t_32>0)
    outcome_32(end-t_32+1:end) = []; %remove the redundent part from outcome vector
end

[~,idx_32] = sort(outcome_32,'descend'); 
Individual_32 = Individual_32(idx_32);
V_result_32(g_32) = Individual_32(1); %record the vaccinated unit
Individual_32(1) = []; %drop the unit with highest marginal gain from ground set
Individual_32 = sort(Individual_32);
V_32(V_result_32(g_32),g_32:end) = 1;
t_32 = t_32+1;
g_32 = g_32+1;
end

%start iteration for N_3, edges = full
while (g_33 < d_3+1)
i = 1;
while (i < N_3-t_33+1)
    V_33(Individual_33(i),g_33) = 1;
    F_33(i,g_33) = obj(V_33(:,g_33),W_9,C_3);
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
V_result_11 = nonzeros(V_result_11);
V_result_12 = nonzeros(V_result_12);
V_result_13 = nonzeros(V_result_13);
V_result_21 = nonzeros(V_result_21);
V_result_22 = nonzeros(V_result_22);
V_result_23 = nonzeros(V_result_23);
V_result_31 = nonzeros(V_result_31);
V_result_32 = nonzeros(V_result_32);
V_result_33 = nonzeros(V_result_33);
%sort the vaccinated unit
V_result_11 = sort(V_result_11);
V_result_12 = sort(V_result_12);
V_result_13 = sort(V_result_13);
V_result_21 = sort(V_result_21);
V_result_22 = sort(V_result_22);
V_result_23 = sort(V_result_23);
V_result_31 = sort(V_result_31);
V_result_32 = sort(V_result_32);
V_result_33 = sort(V_result_33);
%save the final result of social security
F_result_11 = max(F_11(:,d_1));
F_result_12 = max(F_12(:,d_1));
F_result_13 = max(F_13(:,d_1));
F_result_21 = max(F_21(:,d_2));
F_result_22 = max(F_22(:,d_2));
F_result_23 = max(F_23(:,d_2));
F_result_31 = max(F_31(:,d_3));
F_result_32 = max(F_32(:,d_3));
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
const_11 = 0;
while (nn <= N_1)
while (mm <= N_1)
   const_11 = const_11+(beta_11*a_1(nn)*a_1(mm)+beta_12*a_1(nn)*b_1(mm)+beta_21*b_1(nn)*a_1(mm)+ beta_22*b_1(nn)*b_1(mm))*adj_1(nn,mm)*I_1(mm)*S_1(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_1
const_11 = const_11/N_1;

mm = 1;
nn = 1;
const_12 = 0;
while (nn <= N_1)
while (mm <= N_1)
   const_12 = const_12+(beta_11*a_1(nn)*a_1(mm)+beta_12*a_1(nn)*b_1(mm)+beta_21*b_1(nn)*a_1(mm)+ beta_22*b_1(nn)*b_1(mm))*adj_2(nn,mm)*I_1(mm)*S_1(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_2
const_12 = const_12/N_1;

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
const_21 = 0;
while (nn <= N_2)
while (mm <= N_2)
   const_21 = const_21+(beta_11*a_2(nn)*a_2(mm)+beta_12*a_2(nn)*b_2(mm)+beta_21*b_2(nn)*a_2(mm)+ beta_22*b_2(nn)*b_2(mm))*adj_4(nn,mm)*I_2(mm)*S_2(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_4
const_21 = const_21/N_2;

mm = 1;
nn = 1;
const_22 = 0;
while (nn <= N_2)
while (mm <= N_2)
   const_22 = const_22+(beta_11*a_2(nn)*a_2(mm)+beta_12*a_2(nn)*b_2(mm)+beta_21*b_2(nn)*a_2(mm)+ beta_22*b_2(nn)*b_2(mm))*adj_5(nn,mm)*I_2(mm)*S_2(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_5
const_22 = const_22/N_2;

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
const_31 = 0;
while (nn <= N_3)
while (mm <= N_3)
   const_31 = const_31+(beta_11*a_3(nn)*a_3(mm)+beta_12*a_3(nn)*b_3(mm)+beta_21*b_3(nn)*a_3(mm)+ beta_22*b_3(nn)*b_3(mm))*adj_7(nn,mm)*I_3(mm)*S_3(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_7
const_31 = const_31/N_3;

mm = 1;
nn = 1;
const_32 = 0;
while (nn <= N_3)
while (mm <= N_3)
   const_32 = const_32+(beta_11*a_3(nn)*a_3(mm)+beta_12*a_3(nn)*b_3(mm)+beta_21*b_3(nn)*a_3(mm)+ beta_22*b_2(nn)*b_3(mm))*adj_8(nn,mm)*I_3(mm)*S_3(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_8
const_32 = const_32/N_3;

mm = 1;
nn = 1;
const_33 = 0;
while (nn <= N_3)
while (mm <= N_3)
   const_33 = const_33+(beta_11*a_3(nn)*a_3(mm)+beta_12*a_3(nn)*b_3(mm)+beta_21*b_3(nn)*a_3(mm)+ beta_22*b_2(nn)*b_3(mm))*adj_9(nn,mm)*I_3(mm)*S_3(nn);
    mm = mm+1;
end
nn = nn+1;
end
%take average value for interacted constant part with adj_9
const_33 = const_33/N_3;
%combime the previous result together
Y_11 = const_1-const_11+F_result_11;
Y_12 = const_1-const_12+F_result_12;
Y_13 = const_1-const_13+F_result_13;
Y_21 = const_2-const_21+F_result_21;
Y_22 = const_2-const_22+F_result_22;
Y_23 = const_2-const_23+F_result_23;
Y_31 = const_3-const_31+F_result_31;
Y_32 = const_3-const_32+F_result_32;
Y_33 = const_3-const_33+F_result_33;
toc

%clear the irrelevant variables
clear Individual_11 Individual_12 Individual_13 Individual_21 Individual_22 Individual_23 Individual_31 Individual_32 Individual_33
clear outcome_11 outcome_12 outcome_13 outcome_21 outcome_22 outcome_23 outcome_31 outcome_32 outcome_33
clear idx_11 idx_12 idx_13 idx_21 idx_22 idx_23 idx_31 idx_32 idx_33
clear t_11 t_12 t_13 t_21 t_22 t_23 t_31 t_32 t_33
clear i mm nn zz
clear g_11 g_12 g_13 g_21 g_22 g_23 g_31 g_32 g_33
clear V_11 V_12 V_13 V_21 V_22 V_23 V_31 V_32 V_33
clear W_11 W_12 W_13 W_14 W_15 W_16 W_17 W_18 W_19 W_21 W_22 W_23 W_24 W_25 W_26 W_27 W_28 W_29 
clear W_31 W_32 W_33 W_34 W_35 W_36 W_37 W_38 W_39 W_41 W_42 W_43 W_44 W_45 W_46 W_47 W_48 W_49 
clear W_1 W_2 W_3 W_4 W_5 W_6 W_7 W_8 W_9
clear beta_11 beta_12 beta_21 beta_22 gamma_1 gamma_2
clear a_1 a_2 a_3 b_1 b_2 b_3 S_1 S_2 S_3 I_1 I_2 I_3 R_1 R_2 R_3
clear const_1 const_2 const_3 const_11 const_12 const_13 const_21 const_22 const_23 const_31 const_32 const_33
clear F_result_11 F_result_12 F_result_13 F_result_21 F_result_22 F_result_23 F_result_31 F_result_32 F_result_33
clear F_11 F_12 F_13 F_21 F_22 F_23 F_31 F_32 F_33 C_1 C_2 C_3
clear adj_1 adj_2 adj_3 adj_4 adj_5 adj_6 adj_7 adj_8 adj_9
% Save coefficients and other variables for graphs
save(filename_results);

diary off





