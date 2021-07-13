%Comparision Between Greedy and Targeting Without Network Information
%Author: Toru kitagawa and Guanyi Wang

clear all

color = true;%switch: true use color, false use black-white
if (color)
    plot_1_color = [0 0.4470 0.7410];
    plot_2_color = [0.6350 0.0780 0.1840];
    plot_3_color = [0.9290 0.6940 0.1250];
   save_name_500 =  'greedy_nonet_500_color.pdf';%file name for N=500, color
   save_name_800 =  'greedy_nonet_800_color.pdf';%file name for N=800, color
else
    plot_1_color = [0.2 0.2 0.2];
    plot_2_color = [0.5 0.5 0.5];
    plot_3_color = [0.8 0.8 0.8];
   save_name_500 =  'greedy_nonet_500_black.pdf';%file name for N=500, Black
   save_name_800 =  'greedy_nonet_800_black.pdf';%file name for N=800, Black
end


% DATA INPUT AND PREPARATION
load('../greedy_multi_fixden_para1_bgt1.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y21_1_1 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 1 and budget constraint 1, N=500
Y31_1_1 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 1 and budget constraint 1, N=800
load('../greedy_multi_fixden_para2_bgt1.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y22_1_1 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 2 and budget constraint 1, N=500
Y32_1_1 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 2 and budget constraint 1, N=800
load('../nonet_fixden_para1_bgt1.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y21_2_1 = [Y_21_2,Y_22_2,Y_23_2];%save the node for TWNI, parameter 1 and budget constraint 1, N=500
Y31_2_1 = [Y_31_2,Y_32_2,Y_33_2];%save the node for TWNI, parameter 1 and budget constraint 1, N=800
load('../nonet_fixden_para2_bgt1.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y22_2_1 = [Y_21_2,Y_22_2,Y_23_2];%save the node for TWNI, parameter 2 and budget constraint 1, N=500
Y32_2_1 = [Y_31_2,Y_32_2,Y_33_2];%save the node for TWNI, parameter 2 and budget constraint 1, N=800
load('../greedy_multi_fixden_para1_bgt2.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y21_1_2 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 1 and budget constraint 2, N=500
Y31_1_2 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 1 and budget constraint 2, N=800
load('../greedy_multi_fixden_para2_bgt2.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y22_1_2 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 2 and budget constraint 2, N=500
Y32_1_2 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 2 and budget constraint 2, N=800
load('../nonet_fixden_para1_bgt2.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y21_2_2 = [Y_21_2,Y_22_2,Y_23_2];%save the node for TWNI, parameter 1 and budget constraint 2, N=500
Y31_2_2 = [Y_31_2,Y_32_2,Y_33_2];%save the node for TWNI, parameter 1 and budget constraint 2, N=800
load('../nonet_fixden_para2_bgt2.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y22_2_2 = [Y_21_2,Y_22_2,Y_23_2];%save the node for TWNI, parameter 2 and budget constraint 2, N=500
Y32_2_2 = [Y_31_2,Y_32_2,Y_33_2];%save the node for TWNI, parameter 2 and budget constraint 2, N=800
load('../greedy_multi_fixden_para1_bgt3.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y21_1_3 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 1 and budget constraint 3, N=500
Y31_1_3 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 1 and budget constraint 3, N=800
load('../greedy_multi_fixden_para2_bgt3.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y22_1_3 = [Y_21_1,Y_22_1,Y_23_1];%save the node for greedy algorithm, parameter 2 and budget constraint 3, N=500
Y32_1_3 = [Y_31_1,Y_32_1,Y_33_1];%save the node for greedy algorithm, parameter 2 and budget constraint 3, N=800
load('../nonet_fixden_para1_bgt3.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y21_2_3 = [Y_21_2,Y_22_2,Y_23_2];%save the node for TWNI, parameter 1 and budget constraint 3, N=500
Y31_2_3 = [Y_31_2,Y_32_2,Y_33_2];%save the node for TWNI, parameter 1 and budget constraint 3, N=800
load('../nonet_fixden_para2_bgt3.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y22_2_3 = [Y_21_2,Y_22_2,Y_23_2];%save the node for TWNI, parameter 2 and budget constraint 3, N=500
Y32_2_3 = [Y_31_2,Y_32_2,Y_33_2];%save the node for TWMI, parameter 2 and budget constraint 3, N=800
%plot the comparision
sz_1=40; %choose the size of scatter
sz_2=40; %choose the size of scatter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the first one with N=500
figure(1)
%plot the node for density=0.1 (x-axis is the result from greedy algorithm, and y-axis is the result from TWNI)
scatter(Y21_1_1(1), Y21_2_1(1), sz_1,'o', 'MarkerEdgeColor',plot_1_color ,'LineWidth',1.5)%plot the node for N=500
hold on;
scatter(Y21_1_2(1), Y21_2_2(1), sz_1,'o', 'MarkerEdgeColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y21_1_3(1), Y21_2_3(1), sz_1,'o', 'MarkerEdgeColor',plot_3_color ,'LineWidth',1.5)
hold on;
scatter(Y22_1_1(1), Y22_2_1(1), sz_1,'o', 'MarkerEdgeColor',plot_1_color ,'MarkerFaceColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y22_1_2(1), Y22_2_2(1), sz_1,'o', 'MarkerEdgeColor',plot_2_color ,'MarkerFaceColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y22_1_3(1), Y22_2_3(1), sz_1,'o', 'MarkerEdgeColor',plot_3_color ,'MarkerFaceColor',plot_3_color ,'LineWidth',1.5)

hold on;
x_45_1 = linspace(0.2,0.79);%generate nodes for 45 degree line
y_45_1 = linspace(0.2,0.79);%generate nodes for 45 degree line
%draw the 45 degree line
plot(x_45_1,y_45_1,':','color','k');

hold on;
%plot the node for density=1 (x-axis is the result from greedy algorithm, and y-axis is the result from TWNI)
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

% add annotation for density=1
dim_11 = [.6 .58 .28 .25];%the location of rectangle for density=1 in figure 1
dim_12 = [.17 .11 .35 .3];%the location of rectangle for density=1 in figure 2
annotation('rectangle',dim_11,'Color','red')%circle the area of density=1
x_edge_full = [0.5 0.6];
y_edge_full = [0.3 0.55];
annotation('textarrow',x_edge_full,y_edge_full,'String','density = 1');%add word 'edges = full'

% add annotation for density=0.1
dim_21 = [.55 .55 .35 .32];%the location of rectangle for density=0.1
dim_22 = [.6 .6 .3 .32];%the location of rectangle for density=0.1
annotation('rectangle',dim_21,'Color','red')%circle the area of density = 0.1
x_dens_05 = [0.7 0.75];
y_dens_05 = [0.45 0.64];
annotation('textarrow',x_dens_05,y_dens_05,'String','density = 0.1');%add word 'density=0.1'

%add label
legend('N=500,d=7%,  Para=1','N=500,d=10%,Para=1','N=500,d=20%,Para=1','N=500,d=7%,  Para=2','N=500,d=10%,Para=2','N=500,d=20%,Para=2','45-degree line','Location','northwest')
xlabel('Greedy Algorithm')
ylabel('Targeting Without Network Information')
ax = gca;
exportgraphics(ax,save_name_500,'Resolution',300) %save the graph of figure 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the second one with N=800
figure(2)
%plot the node for density=0.1 (x-axis is the result from greedy algorithm, and y-axis is the result from TWNI)
scatter(Y31_1_1(1), Y31_2_1(1), sz_2,'o', 'MarkerEdgeColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y31_1_2(1), Y31_2_2(1), sz_2,'o', 'MarkerEdgeColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y31_1_3(1), Y31_2_3(1), sz_2,'o', 'MarkerEdgeColor',plot_3_color ,'LineWidth',1.5)
hold on;
scatter(Y32_1_1(1), Y32_2_1(1), sz_2,'o', 'MarkerEdgeColor',plot_1_color ,'MarkerFaceColor',plot_1_color ,'LineWidth',1.5)
hold on;
scatter(Y32_1_2(1), Y32_2_2(1), sz_2,'o', 'MarkerEdgeColor',plot_2_color ,'MarkerFaceColor',plot_2_color ,'LineWidth',1.5)
hold on;
scatter(Y32_1_3(1), Y32_2_3(1), sz_2,'o', 'MarkerEdgeColor',plot_3_color ,'MarkerFaceColor',plot_3_color ,'LineWidth',1.5)

hold on;
x_45_2 = linspace(0.2,0.79);
y_45_2 = linspace(0.2,0.79);
%draw the 45 degree line
plot(x_45_2,y_45_2,':','color','k');

hold on;
%plot the node for density=1 (x-axis is the result from greedy algorithm, and y-axis is the result from TWNI)
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

% add annotation for density=1
annotation('rectangle',dim_12,'Color','red')%circle the area of density=1
x_edge_full = [0.55 0.45];
y_edge_full = [0.2 0.3];
annotation('textarrow',x_edge_full,y_edge_full,'String','density = 1');%add word 'edges = full'

% add annotation for density=0.1
annotation('rectangle',dim_21,'Color','red')%circle the area of density=0.1
x_edge_10N = [0.7 Y21_1_1(1)+0.1];
y_edge_10N = [0.4 Y21_2_1(1)];
annotation('textarrow',x_edge_10N,y_edge_10N,'String','density = 0.1');%add word 'edges = 10N'

%add label
legend('N=800,d=7%,  Para=1','N=800,d=10%,Para=1','N=800,d=20%,Para=1','N=800,d=7%,  Para=2','N=800,d=10%,Para=2','N=800,d=20%,Para=2','45-degree line','Location','northwest')
xlabel('Greedy Algorithm')
ylabel('Targeting Without Network Information')
bx = gca;
exportgraphics(bx,save_name_800,'Resolution',300)%save the graph of figure 2