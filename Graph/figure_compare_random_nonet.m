%Comparision Between Random Allocation and Targeting Without Network Information
%Author: Toru kitagawa and Guanyi Wang

clear all

color = false;%switch: true use color, false use black-white
if (color)
    plot_1_color = [0 0.4470 0.7410];
    plot_2_color = [0.6350 0.0780 0.1840];
    plot_3_color = [0.9290 0.6940 0.1250];
   save_name_500 =  'random_nonet_500_color.pdf';
   save_name_800 =  'random_nonet_800_color.pdf';
else
    plot_1_color = [0.2 0.2 0.2];
    plot_2_color = [0.5 0.5 0.5];
    plot_3_color = [0.8 0.8 0.8];
   save_name_500 =  'random_nonet_500_black.pdf';
   save_name_800 =  'random_nonet_800_black.pdf';
end
% DATA INPUT AND PREPARATION
load('../random_multi_para1_bgt1.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y21_1_1 = [Y_21_1,Y_22_1,Y_23_1];
Y31_1_1 = [Y_31_1,Y_32_1,Y_33_1];
load('../random_multi_para2_bgt1.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y22_1_1 = [Y_21_1,Y_22_1,Y_23_1];
Y32_1_1 = [Y_31_1,Y_32_1,Y_33_1];
load('../nonet_multi_para1_bgt1.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y21_2_1 = [Y_21_2,Y_22_2,Y_23_2];
Y31_2_1 = [Y_31_2,Y_32_2,Y_33_2];
load('../nonet_multi_para2_bgt1.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y22_2_1 = [Y_21_2,Y_22_2,Y_23_2];
Y32_2_1 = [Y_31_2,Y_32_2,Y_33_2];
load('../random_multi_para1_bgt2.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y21_1_2 = [Y_21_1,Y_22_1,Y_23_1];
Y31_1_2 = [Y_31_1,Y_32_1,Y_33_1];
load('../random_multi_para2_bgt2.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y22_1_2 = [Y_21_1,Y_22_1,Y_23_1];
Y32_1_2 = [Y_31_1,Y_32_1,Y_33_1];
load('../nonet_multi_para1_bgt2.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y21_2_2 = [Y_21_2,Y_22_2,Y_23_2];
Y31_2_2 = [Y_31_2,Y_32_2,Y_33_2];
load('../nonet_multi_para2_bgt2.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y22_2_2 = [Y_21_2,Y_22_2,Y_23_2];
Y32_2_2 = [Y_31_2,Y_32_2,Y_33_2];
load('../random_multi_para1_bgt3.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y21_1_3 = [Y_21_1,Y_22_1,Y_23_1];
Y31_1_3 = [Y_31_1,Y_32_1,Y_33_1];
load('../random_multi_para2_bgt3.mat');
Y_21_1 = Y_21_ave;
Y_22_1 = Y_22_ave;
Y_23_1 = Y_23_ave;
Y_31_1 = Y_31_ave;
Y_32_1 = Y_32_ave;
Y_33_1 = Y_33_ave;
Y22_1_3 = [Y_21_1,Y_22_1,Y_23_1];
Y32_1_3 = [Y_31_1,Y_32_1,Y_33_1];
load('../nonet_multi_para1_bgt3.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y21_2_3 = [Y_21_2,Y_22_2,Y_23_2];
Y31_2_3 = [Y_31_2,Y_32_2,Y_33_2];
load('../nonet_multi_para2_bgt3.mat');
Y_21_2 = Y_21_ave;
Y_22_2 = Y_22_ave;
Y_23_2 = Y_23_ave;
Y_31_2 = Y_31_ave;
Y_32_2 = Y_32_ave;
Y_33_2 = Y_33_ave;
Y22_2_3 = [Y_21_2,Y_22_2,Y_23_2];
Y32_2_3 = [Y_31_2,Y_32_2,Y_33_2];
%plot the comparision
sz_1=40; %choose the size of scatter
sz_2=40; %choose the size of scatter
%plot the first one with N=500
figure(1)
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
x_45_1 = linspace(0,0.8);
y_45_1 = linspace(0,0.8);
plot(x_45_1,y_45_1,':','color','k');
hold on;
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
% add annotation
dim_11 = [.3 .28 .33 .31];
dim_12 = [.13 .11 .35 .3];
annotation('rectangle',dim_11,'Color','red')
x_edge_full = [0.7 Y21_1_1(3)+0.2];
y_edge_full = [0.2 Y21_2_1(3)];
annotation('textarrow',x_edge_full,y_edge_full,'String','edges = full');
% add annotation
dim_2 = [.6 .6 .3 .32];
annotation('rectangle',dim_2,'Color','red')
x_edge_10N = [0.7 Y21_1_1(2)+0.1];
y_edge_10N = [0.4 Y21_2_1(2)];
annotation('textarrow',x_edge_10N,y_edge_10N,'String','edges = 10N');
%add label
legend('N=500,d=7%,  Para=1','N=500,d=10%,Para=1','N=500,d=20%,Para=1','N=500,d=7%,  Para=2','N=500,d=10%,Para=2','N=500,d=20%,Para=2','45-degree line','Location','northwest')
xlabel('Random Algorithm')
ylabel('Targeting Without Network Information')
ax = gca;
exportgraphics(ax,save_name_500,'Resolution',300) 
%plot the second one with N=800
figure(2)
%plot the second one with N=800
figure(2)
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
x_45_2 = linspace(0,0.8);
y_45_2 = linspace(0,0.8);
plot(x_45_2,y_45_2,':','color','k');
hold on;
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
% add annotation
annotation('rectangle',dim_12,'Color','red')
x_edge_full = [0.7 Y21_1_1(3)+0.2];
y_edge_full = [0.2 Y21_2_1(3)];
annotation('textarrow',x_edge_full,y_edge_full,'String','edges = full');
% add annotation
annotation('rectangle',dim_2,'Color','red')
x_edge_10N = [0.7 Y21_1_1(2)+0.1];
y_edge_10N = [0.4 Y21_2_1(2)];
annotation('textarrow',x_edge_10N,y_edge_10N,'String','edges = 10N');
%add label
legend('N=800,d=7%,  Para=1','N=800,d=10%,Para=1','N=800,d=20%,Para=1','N=800,d=7%,  Para=2','N=800,d=10%,Para=2','N=800,d=20%,Para=2','45-degree line','Location','northwest')
xlabel('Random Algorithm')
ylabel('Targeting Without Network Information')
bx = gca;
exportgraphics(bx,save_name_800,'Resolution',300)