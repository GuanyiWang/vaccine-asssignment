%Brute-force search the global optimization
%Author: Guanyi Wang and Toru Kitagawa

clear all

% DATA INPUT AND PREPARATION
load('../simulation_small.mat')

N_1 = length(S_1); 
N_2 = length(S_2);
N_3 = length(S_3);

% Parameter setting 
parameter = false; % switch: choose 1st parameter set into account or 2nd one
budget_constrant = 3; %swtich: choose the capacity constraint level from 1, 2, 3

if (budget_constrant == 1)
if (parameter)
diary('global_para1.log');
filename_results = '../global_para1_bgt1.mat';
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
diary('global_multi_para2.log');
filename_results = '../global_para2_bgt1.mat';
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
if (budget_constrant == 2)
    if (parameter)
diary('global_para1.log');
filename_results = '../global_para1_bgt2.mat';
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
diary('global_para2.log');
filename_results = '../global_para2_bgt2.mat';
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

if (budget_constrant == 3)
    if (parameter)
diary('global_para1.log');
filename_results = '../global_para1_bgt3.mat';
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
diary('global_para2.log');
filename_results = '../global_para2_bgt3.mat';
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


%Define the parameter matrix
C_1 = ones(N_1,1) - R_1 - diag(a_1*I_1')*gamma_1 - diag(b_1*I_1')*gamma_2 - S_1; %define the vector C_1
C_2 = ones(N_2,1) - R_2 - diag(a_2*I_2')*gamma_1 - diag(b_2*I_2')*gamma_2 - S_2; %define the vector C_2
C_3 = ones(N_3,1) - R_3 - diag(a_3*I_3')*gamma_1 - diag(b_3*I_3')*gamma_2 - S_3; %define the vector C_3


W_11 = beta_11*diag(a_1)*diag(S_1)*adj_1*diag(a_1)*diag(I_1);
W_12 = beta_11*diag(a_1)*diag(S_1)*adj_2*diag(a_1)*diag(I_1);
W_13 = beta_11*diag(a_1)*diag(S_1)*adj_3*diag(a_1)*diag(I_1);
W_14 = beta_11*diag(a_2)*diag(S_2)*adj_4*diag(a_2)*diag(I_2);
W_15 = beta_11*diag(a_2)*diag(S_2)*adj_5*diag(a_2)*diag(I_2);
W_16 = beta_11*diag(a_2)*diag(S_2)*adj_6*diag(a_2)*diag(I_2);
W_17 = beta_11*diag(a_3)*diag(S_3)*adj_7*diag(a_3)*diag(I_3);
W_18 = beta_11*diag(a_3)*diag(S_3)*adj_8*diag(a_3)*diag(I_3);
W_19 = beta_11*diag(a_3)*diag(S_3)*adj_9*diag(a_3)*diag(I_3);

W_21 = beta_12*diag(a_1)*diag(S_1)*adj_1*diag(b_1)*diag(I_1);
W_22 = beta_12*diag(a_1)*diag(S_1)*adj_2*diag(b_1)*diag(I_1);
W_23 = beta_12*diag(a_1)*diag(S_1)*adj_3*diag(b_1)*diag(I_1);
W_24 = beta_12*diag(a_2)*diag(S_2)*adj_4*diag(b_2)*diag(I_2);
W_25 = beta_12*diag(a_2)*diag(S_2)*adj_5*diag(b_2)*diag(I_2);
W_26 = beta_12*diag(a_2)*diag(S_2)*adj_6*diag(b_2)*diag(I_2);
W_27 = beta_12*diag(a_3)*diag(S_3)*adj_7*diag(b_3)*diag(I_3);
W_28 = beta_12*diag(a_3)*diag(S_3)*adj_8*diag(b_3)*diag(I_3);
W_29 = beta_12*diag(a_3)*diag(S_3)*adj_9*diag(b_3)*diag(I_3);

W_31 = beta_21*diag(b_1)*diag(S_1)*adj_1*diag(a_1)*diag(I_1);
W_32 = beta_21*diag(b_1)*diag(S_1)*adj_2*diag(a_1)*diag(I_1);
W_33 = beta_21*diag(b_1)*diag(S_1)*adj_3*diag(a_1)*diag(I_1);
W_34 = beta_21*diag(b_2)*diag(S_2)*adj_4*diag(a_2)*diag(I_2);
W_35 = beta_21*diag(b_2)*diag(S_2)*adj_5*diag(a_2)*diag(I_2);
W_36 = beta_21*diag(b_2)*diag(S_2)*adj_6*diag(a_2)*diag(I_2);
W_37 = beta_21*diag(b_3)*diag(S_3)*adj_7*diag(a_3)*diag(I_3);
W_38 = beta_21*diag(b_3)*diag(S_3)*adj_8*diag(a_3)*diag(I_3);
W_39 = beta_21*diag(b_3)*diag(S_3)*adj_9*diag(a_3)*diag(I_3);

W_41 = beta_22*diag(b_1)*diag(S_1)*adj_1*diag(b_1)*diag(I_1);
W_42 = beta_22*diag(b_1)*diag(S_1)*adj_2*diag(b_1)*diag(I_1);
W_43 = beta_22*diag(b_1)*diag(S_1)*adj_3*diag(b_1)*diag(I_1);
W_44 = beta_22*diag(b_2)*diag(S_2)*adj_4*diag(b_2)*diag(I_2);
W_45 = beta_22*diag(b_2)*diag(S_2)*adj_5*diag(b_2)*diag(I_2);
W_46 = beta_22*diag(b_2)*diag(S_2)*adj_6*diag(b_2)*diag(I_2);
W_47 = beta_22*diag(b_3)*diag(S_3)*adj_7*diag(b_3)*diag(I_3);
W_48 = beta_22*diag(b_3)*diag(S_3)*adj_8*diag(b_3)*diag(I_3);
W_49 = beta_22*diag(b_3)*diag(S_3)*adj_9*diag(b_3)*diag(I_3);
%define the weighting matrix W
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
%brute-force search
Individual_1 = (1:N_1);
Individual_2 = (1:N_2);
Individual_3 = (1:N_3);

%find all the possible combination of vaccine choice in N_1
Ind_1 = nchoosek(Individual_1,d_1); 
Ind_1 = Ind_1';
%find all the possible combination of vaccine choice in N_2
Ind_2 = nchoosek(Individual_2,d_2);
Ind_2 = Ind_2';
%find all the possible combination of vaccine choice in N_3
Ind_3 = nchoosek(Individual_3,d_3); 
Ind_3 = Ind_3';
%find the total number of possible combinations
eend_1 = nchoosek(N_1,d_1);
eend_2 = nchoosek(N_2,d_2);
eend_3 = nchoosek(N_3,d_3);
%Initialize the vector of V
V_11 = zeros(N_1,eend_1); 
V_12 = zeros(N_1,eend_1); 
V_13 = zeros(N_1,eend_1); 
V_21 = zeros(N_2,eend_2); 
V_22 = zeros(N_2,eend_2); 
V_23 = zeros(N_2,eend_2);
V_31 = zeros(N_3,eend_3); 
V_32 = zeros(N_3,eend_3); 
V_33 = zeros(N_3,eend_3);
%Set the vector of objective function
F_11 = zeros(eend_1,1); 
F_12 = zeros(eend_1,1); 
F_13 = zeros(eend_1,1); 
F_21 = zeros(eend_2,1);
F_22 = zeros(eend_2,1);
F_23 = zeros(eend_2,1);
F_31 = zeros(eend_3,1);
F_32 = zeros(eend_3,1);
F_33 = zeros(eend_3,1);
tic
ii = 1; %initialize the paramter
while (ii<=eend_1)
    i= 1;
while (i<=d_1)
V_11(Ind_1(i,ii),ii) = 1;
i = i+1;
end
    F_11(ii) = obj(V_11(:,ii),W_1,C_1);
    ii = ii+1;
end

ii = 1; %initialize the paramter
while (ii<=eend_1)
    i= 1;
while (i<=d_1)
V_12(Ind_1(i,ii),ii) = 1;
i = i+1;
end
    F_12(ii) = obj(V_12(:,ii),W_2,C_1);
    ii = ii+1;
end

ii = 1; %initialize the paramter
while (ii<=eend_1)
    i= 1;
while (i<=d_1)
V_13(Ind_1(i,ii),ii) = 1;
i = i+1;
end
    F_13(ii) = obj(V_13(:,ii),W_3,C_1);
    ii = ii+1;
end

ii = 1; %initialize the paramter
while (ii<=eend_2)
    i= 1;
while (i<=d_2)
V_21(Ind_2(i,ii),ii) = 1;
i = i+1;
end
    F_21(ii) = obj(V_21(:,ii),W_4,C_2);
    ii = ii+1;
end

ii = 1; %initialize the paramter
while (ii<=eend_2)
    i= 1;
while (i<=d_2)
V_22(Ind_2(i,ii),ii) = 1;
i = i+1;
end
    F_22(ii) = obj(V_22(:,ii),W_5,C_2);
    ii = ii+1;
end

ii = 1; %initialize the paramter
while (ii<=eend_2)
    i= 1;
while (i<=d_2)
V_23(Ind_2(i,ii),ii) = 1;
i = i+1;
end
    F_23(ii) = obj(V_23(:,ii),W_6,C_2);
    ii = ii+1;
end

ii = 1; %initialize the paramter
while (ii<=eend_3)
    i= 1;
while (i<=d_3)
V_31(Ind_3(i,ii),ii) = 1;
i = i+1;
end
    F_31(ii) = obj(V_31(:,ii),W_7,C_3);
    ii = ii+1;
end

ii = 1; %initialize the paramter
while (ii<=eend_3)
    i= 1;
while (i<=d_3)
V_32(Ind_3(i,ii),ii) = 1;
i = i+1;
end
    F_32(ii) = obj(V_32(:,ii),W_8,C_3);
    ii = ii+1;
end

ii = 1; %initialize the paramter
while (ii<=eend_3)
    i= 1;
while (i<=d_3)
V_33(Ind_3(i,ii),ii) = 1;
i = i+1;
end
    F_33(ii) = obj(V_33(:,ii),W_9,C_3);
    ii = ii+1;
end
[F_result_11,idx_11] = max(F_11);
[F_result_12,idx_12] = max(F_12);
[F_result_13,idx_13] = max(F_13);
[F_result_21,idx_21] = max(F_21);
[F_result_22,idx_22] = max(F_22);
[F_result_23,idx_23] = max(F_23);
[F_result_31,idx_31] = max(F_31);
[F_result_32,idx_32] = max(F_32);
[F_result_33,idx_33] = max(F_33);

%return the index of vaccinated agents
V_output_11 = Ind_1(:,idx_11);
V_output_12 = Ind_1(:,idx_12);
V_output_13 = Ind_1(:,idx_13);
V_output_21 = Ind_2(:,idx_21);
V_output_22 = Ind_2(:,idx_22);
V_output_23 = Ind_2(:,idx_23);
V_output_31 = Ind_3(:,idx_31);
V_output_32 = Ind_3(:,idx_32);
V_output_33 = Ind_3(:,idx_33);
%calculate the initial objective function value
zz = 1;
const_1 = 0;
while (zz<=N_1)
   const_1 = const_1 + R_1(zz)+a_1(zz)*I_1(zz)*gamma_1+b_1(zz)*I_1(zz)*gamma_2+S_1(zz);
    zz = zz+1;
end
%find the average of first part
const_1 = const_1/N_1;
%calcuate the interacted constant part sepearately 
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

const_13 = const_13/N_1;

zz = 1;
const_2 = 0;
while (zz<=N_2)
   const_2 = const_2 + R_2(zz)+a_2(zz)*I_2(zz)*gamma_2+b_2(zz)*I_2(zz)*gamma_2+S_2(zz);
    zz = zz+1;
end

const_2 = const_2/N_2;

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

const_23 = const_23/N_2;

zz = 1;
const_3 = 0;
while (zz<=N_3)
   const_3 = const_3 + R_3(zz)+a_3(zz)*I_3(zz)*gamma_1+b_3(zz)*I_3(zz)*gamma_2+S_3(zz);
    zz = zz+1;
end


const_3 = const_3/N_3;

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
%generate the graph
G_1 = graph(adj_1);
G_2 = graph(adj_2);
G_3 = graph(adj_3);
G_4 = graph(adj_4);
G_5 = graph(adj_5);
G_6 = graph(adj_6);
G_7 = graph(adj_7);
G_8 = graph(adj_8);
G_9 = graph(adj_9);
%plot the specific graph with highlight nodes
%h_4 = plot(G_4);
%highlight(h_4,[V_output_21(1:end)],'NodeColor','r','MarkerSize',6);
%h_4.NodeLabel = {};
%ax = gca;
%exportgraphics(ax,'Global_30.pdf','Resolution',300) 
%h_8 = plot(G_8);
%highlight(h_8,[V_output_32(1:end)],'NodeColor','r','MarkerSize',6);
%h_8.NodeLabel = {};
%bx = gca;
%exportgraphics(bx,'Global_35.pdf','Resolution',300) 
clear Ind_1 Ind_2 Ind_3 eend_1 eend_2 eend_3 C_1 C_2 C_3
clear W_1 W_2 W_3 W_4 W_5 W_6 W_7 W_8 W_9 idx_11 idx_12 idx_13 idx_21 idx_22 idx_23 idx_31 idx_32 idx_33
clear beta_11 beta_12 beta_21 beta_22 gamma_1 gamma_2
clear a_1 a_2 a_3 b_1 b_2 b_3 S_1 S_2 S_3 I_1 I_2 I_3 R_1 R_2 R_3
clear V_11 V_12 V_13 V_21 V_22 V_23 V_31 V_32 V_33
clear W_11 W_12 W_13 W_14 W_15 W_16 W_17 W_18 W_19 W_21 W_22 W_23 W_24 W_25 W_26 W_27 W_28 W_29 
clear W_31 W_32 W_33 W_34 W_35 W_36 W_37 W_38 W_39 W_41 W_42 W_43 W_44 W_45 W_46 W_47 W_48 W_49
clear F_11 F_12 F_13 F_21 F_22 F_23 F_31 F_32 F_33
clear F_result_11 F_result_12 F_result_13 F_result_21 F_result_22 F_result_23 F_result_31 F_result_32 F_result_33
clear i mm nn zz ii cc
clear const_1 const_2 const_3 const_11 const_12 const_13 const_21 const_22 const_23 const_31 const_32 const_33
clear G_1 G_2 G_3 G_4 G_5 G_6 G_7 G_8 G_9 adj_1 adj_2 adj_3 adj_4 adj_5 adj_6 adj_7 adj_8 adj_9
clear Individual_1 Individual_2 Individual_3
% Save coefficients and other variables for graphs
save(filename_results);

diary off