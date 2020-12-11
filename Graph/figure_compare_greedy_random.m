%Comparision Between Greedy and Random Allocation
%Author: Toru kitagawa and Guanyi Wang


clear all

color = false;%switch: true use color, false use black-white
if (color)
    plot_1_color = [0 0.4470 0.7410];
    plot_2_color = [0.6350 0.0780 0.1840];
    plot_3_color = [0.9290 0.6940 0.1250];
   save_name_500 =  'greedy_random_500_color.pdf';%file name for N=500, color
   save_name_800 =  'greedy_random_800_color.pdf';%file name for N=800, color
else
    plot_1_color = [0.2 0.2 0.2];
    plot_2_color = [0.5 0.5 0.5];
    plot_3_color = [0.8 0.8 0.8];
   save_name_500 =  'greedy_random_500_black.pdf';%file name for N=500, Black
   save_name_800 =  'greedy_random_800_black.pdf';%file name for N=800, Black
end
% DATA INPUT AND PREPARATION
load('../greedy_multi_para1_bgt1.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y21_1_1 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 1 and budget constraint 1, N=500
Y31_1_1 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 1 and budget constraint 1, N=800
load('../greedy_multi_para2_bgt1.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y22_1_1 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 2 and budget constraint 1, N=500
Y32_1_1 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 2 and budget constraint 1, N=800
load('../random_multi_para1_bgt1.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y21_2_1 = [Y_21_2,Y_22_2,Y_23_2];%save the node for random allocation, parameter 1 and budget constraint 1, N=500
Y31_2_1 = [Y_31_2,Y_32_2,Y_33_2];%save the node for random allocation, parameter 1 and budget constraint 1, N=800
load('../random_multi_para2_bgt1.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y22_2_1 = [Y_21_2,Y_22_2,Y_23_2];%save the node for random allocation, parameter 2 and budget constraint 1, N=500
Y32_2_1 = [Y_31_2,Y_32_2,Y_33_2];%save the node for random allocation, parameter 2 and budget constraint 1, N=800
load('../greedy_multi_para1_bgt2.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y21_1_2 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 1 and budget constraint 2, N=500
Y31_1_2 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 1 and budget constraint 2, N=800
load('../greedy_multi_para2_bgt2.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y22_1_2 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 2 and budget constraint 2, N=500
Y32_1_2 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 2 and budget constraint 2, N=800
load('../random_multi_para1_bgt2.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y21_2_2 = [Y_21_2,Y_22_2,Y_23_2];%save the node for random allocation, parameter 1 and budget constraint 2, N=500
Y31_2_2 = [Y_31_2,Y_32_2,Y_33_2];%save the node for random allocation, parameter 1 and budget constraint 2, N=800
load('../random_multi_para2_bgt2.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y22_2_2 = [Y_21_2,Y_22_2,Y_23_2];%save the node for random allocation, parameter 2 and budget constraint 2, N=500
Y32_2_2 = [Y_31_2,Y_32_2,Y_33_2];%save the node for random allocation, parameter 2 and budget constraint 2, N=800
load('../greedy_multi_para1_bgt3.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y21_1_3 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 1 and budget constraint 3, N=500
Y31_1_3 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 1 and budget constraint 3, N=800
load('../greedy_multi_para2_bgt3.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y22_1_3 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 2 and budget constraint 3, N=500
Y32_1_3 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 2 and budget constraint 3, N=800
load('../random_multi_para1_bgt3.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y21_2_3 = [Y_21_2,Y_22_2,Y_23_2];%save the node for random allocation, parameter 1 and budget constraint 3, N=500
Y31_2_3 = [Y_31_2,Y_32_2,Y_33_2];%save the node for random allocation, parameter 1 and budget constraint 3, N=800
load('../random_multi_para2_bgt3.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y22_2_3 = [Y_21_2,Y_22_2,Y_23_2];%save the node for random allocation, parameter 2 and budget constraint 3, N=500
Y32_2_3 = [Y_31_2,Y_32_2,Y_33_2];%save the node for random allocation, parameter 2 and budget constraint 3, N=800
%plot the comparision
sz_1=40; %choose the size of scatter of figure 1
sz_2=40; %choose the size of scatter of figure 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the first one with N=500
figure(1)
%plot the node for edges = 10N (x-axis is the result from greedy algorithm, and y-axis is the result from random allocation)
scatter(Y21_1_1(2), Y21_2_1(2), sz_1,'o', 'MarkerEdgeColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y21_1_2(2), Y21_2_2(2), sz_1,'o', 'MarkerEdgeColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y21_1_3(2), Y21_2_3(2), sz_1,'o', 'MarkerEdgeColor',plot_3_color ,'LineWidth',1.5)
hold on;
scatter(Y22_1_1(2), Y22_2_1(2), sz_1,'o', 'MarkerEdgeColor',plot_1_color ,'MarkerFaceColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y22_1_2(2), Y22_2_2(2), sz_1,'o', 'MarkerEdgeColor',plot_2_color ,'MarkerFaceColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y22_1_3(2), Y22_2_3(2), sz_1,'o', 'MarkerEdgeColor',plot_3_color ,'MarkerFaceColor',plot_3_color ,'LineWidth',1.5)

hold on;
x_45_1 = linspace(0,0.8);%generate nodes for 45 degree line
y_45_1 = linspace(0,0.8);%generate nodes for 45 degree line
%draw the 45 degree line
plot(x_45_1,y_45_1,':','color','k');

hold on;
%plot the node for edges = full (x-axis is the result from greedy algorithm, and y-axis is the result from random allocation)
scatter(Y21_1_1(3), Y21_2_1(3), sz_1,'o', 'MarkerEdgeColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y21_1_2(3), Y21_2_2(3), sz_1,'o', 'MarkerEdgeColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y21_1_3(3), Y21_2_3(3), sz_1,'o', 'MarkerEdgeColor',plot_3_color ,'LineWidth',1.5)
hold on;
scatter(Y22_1_1(3), Y22_2_1(3), sz_1,'o', 'MarkerEdgeColor',plot_1_color ,'MarkerFaceColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y22_1_2(3), Y22_2_2(3), sz_1,'o', 'MarkerEdgeColor',plot_2_color ,'MarkerFaceColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y22_1_3(3), Y22_2_3(3), sz_1,'o', 'MarkerEdgeColor',plot_3_color ,'MarkerFaceColor',plot_3_color ,'LineWidth',1.5)

% add annotation for edges = full
dim_11 = [.3 .28 .33 .31];%the location of rectangle for edges = full in figure 1
dim_12 = [.13 .11 .35 .3];%the location of rectangle for edges = full in figure 2
annotation('rectangle',dim_11,'Color','red')%circle the area of edges = full
x_edge_full = [0.7 Y21_1_1(3)+0.2];
y_edge_full = [0.2 Y21_2_1(3)];
annotation('textarrow',x_edge_full,y_edge_full,'String','edges = full');%add word 'edges = full'

% add annotation for edges = 10N
dim_2 = [.6 .6 .3 .32];%the location of rectangle for edges = 10N
annotation('rectangle',dim_2,'Color','red')%circle the area of edges = 10N
x_edge_10N = [0.7 Y21_1_1(2)+0.1];
y_edge_10N = [0.4 Y21_2_1(2)];
annotation('textarrow',x_edge_10N,y_edge_10N,'String','edges = 10N');%add word 'edges = 10N'

%add label
legend('N=500,d=7%,  Para=1','N=500,d=10%,Para=1','N=500,d=20%,Para=1','N=500,d=7%,  Para=2','N=500,d=10%,Para=2','N=500,d=20%,Para=2','45-degree line','Location','northwest')
xlabel('Greedy Algorithm')
ylabel('Random Allocation')
ax = gca;
exportgraphics(ax,save_name_500,'Resolution',300) %save the graph of figure 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the second one with N=800
figure(2)
%plot the node for edges = 10N (x-axis is the result from greedy algorithm, and y-axis is the result from random allocation)
scatter(Y31_1_1(2), Y31_2_1(2), sz_2,'o', 'MarkerEdgeColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y31_1_2(2), Y31_2_2(2), sz_2,'o', 'MarkerEdgeColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y31_1_3(2), Y31_2_3(2), sz_2,'o', 'MarkerEdgeColor',plot_3_color ,'LineWidth',1.5)
hold on;
scatter(Y32_1_1(2), Y32_2_1(2), sz_2,'o', 'MarkerEdgeColor',plot_1_color ,'MarkerFaceColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y32_1_2(2), Y32_2_2(2), sz_2,'o', 'MarkerEdgeColor',plot_2_color ,'MarkerFaceColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y32_1_3(2), Y32_2_3(2), sz_2,'o', 'MarkerEdgeColor',plot_3_color ,'MarkerFaceColor',plot_3_color ,'LineWidth',1.5)

hold on;
x_45_2 = linspace(0,0.8);%generate nodes for 45 degree line
y_45_2 = linspace(0,0.8);%generate nodes for 45 degree line
%draw the 45 degree line
plot(x_45_2,y_45_2,':','color','k');

hold on;
%plot the node for edges = full (x-axis is the result from greedy algorithm, and y-axis is the result from random allocation)
scatter(Y31_1_1(3), Y31_2_1(3), sz_2,'o', 'MarkerEdgeColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y31_1_2(3), Y31_2_2(3), sz_2,'o', 'MarkerEdgeColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y31_1_3(3), Y31_2_3(3), sz_2,'o', 'MarkerEdgeColor',plot_3_color ,'LineWidth',1.5)
hold on;
scatter(Y32_1_1(3), Y32_2_1(3), sz_2,'o', 'MarkerEdgeColor',plot_1_color ,'MarkerFaceColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y32_1_2(3), Y32_2_2(3), sz_2,'o', 'MarkerEdgeColor',plot_2_color ,'MarkerFaceColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y32_1_3(3), Y32_2_3(3), sz_2,'o', 'MarkerEdgeColor',plot_3_color ,'MarkerFaceColor',plot_3_color ,'LineWidth',1.5)

% add annotation for edges = full
annotation('rectangle',dim_12,'Color','red')%circle the area of edges = full
x_edge_full = [0.7 Y21_1_1(3)+0.2];
y_edge_full = [0.2 Y21_2_1(3)];
annotation('textarrow',x_edge_full,y_edge_full,'String','edges = full');%add word 'edges = full'

% add annotation for edges = 10N
annotation('rectangle',dim_2,'Color','red')%circle the area of edges = 10N
x_edge_10N = [0.7 Y21_1_1(2)+0.1];
y_edge_10N = [0.4 Y21_2_1(2)];
annotation('textarrow',x_edge_10N,y_edge_10N,'String','edges = 10N');%add word 'edges = 10N'

%add label
legend('N=800,d=7%,  Para=1','N=800,d=10%,Para=1','N=800,d=20%,Para=1','N=800,d=7%,  Para=2','N=800,d=10%,Para=2','N=800,d=20%,Para=2','45-degree line','Location','northwest')
xlabel('Greedy Algorithm')
ylabel('Random Allocation')
bx = gca;
exportgraphics(bx,save_name_800,'Resolution',300) %save the graph of figure 2