clear all;
close all;


color1=[0, 0.4470, 0.7410];
color2=[0.8500, 0.3250, 0.0980];
color3=[0.9290, 0.6940, 0.1250];
color4=[0.4940, 0.1840, 0.5560];
color5=[0.4660, 0.6740, 0.1880];
color6=[0.3010, 0.7450, 0.9330];
color7=[0.6350, 0.0780, 0.1840];
color8=[0, 0, 0];

% Import Nup62c run1
Nup62c11=importdata('Nup62_Nup58_Nup54_32_surf1.txt');
Nup62c11=Nup62c11.data;
Nup62c12=importdata('Nup62_Nup58_Nup54_32_surf2.txt');
Nup62c12=Nup62c12.data;
% Import Nup62c run2
Nup62c21=importdata('Nup62_Nup58_Nup54_32_2_surf1.txt');
Nup62c21=Nup62c21.data;
Nup62c22=importdata('Nup62_Nup58_Nup54_32_2_surf2.txt');
Nup62c22=Nup62c22.data;
% Import Nup62c run3
Nup62c31=importdata('Nup62_Nup58_Nup54_32_3_surf1.txt');
Nup62c31=Nup62c31.data;
Nup62c32=importdata('Nup62_Nup58_Nup54_32_3_surf2.txt');
Nup62c32=Nup62c32.data;
% Import Nup62c run4
Nup62c41=importdata('Nup62_Nup58_Nup54_32_4_surf1.txt');
Nup62c41=Nup62c41.data;
Nup62c42=importdata('Nup62_Nup58_Nup54_32_4_surf2.txt');
Nup62c42=Nup62c42.data;
% Import Nup62c run5
Nup62c51=importdata('Nup62_Nup58_Nup54_32_5_surf1.txt');
Nup62c51=Nup62c51.data;
Nup62c52=importdata('Nup62_Nup58_Nup54_32_5_surf2.txt');
Nup62c52=Nup62c52.data;
% Import Nup62c run6
Nup62c61=importdata('Nup62_Nup58_Nup54_32_6_surf1.txt');
Nup62c61=Nup62c61.data;
Nup62c62=importdata('Nup62_Nup58_Nup54_32_6_surf2.txt');
Nup62c62=Nup62c62.data;
% Import Nup62c run7
Nup62c71=importdata('Nup62_Nup58_Nup54_32_7_surf1.txt');
Nup62c71=Nup62c71.data;
Nup62c72=importdata('Nup62_Nup58_Nup54_32_7_surf2.txt');
Nup62c72=Nup62c72.data;
% Import Nup62c run8
Nup62c81=importdata('Nup62_Nup58_Nup54_32_8_surf1.txt');
Nup62c81=Nup62c81.data;
Nup62c82=importdata('Nup62_Nup58_Nup54_32_8_surf2.txt');
Nup62c82=Nup62c82.data;
% Import Nup62c run9
Nup62c91=importdata('Nup62_Nup58_Nup54_32_9_surf1.txt');
Nup62c91=Nup62c91.data;
Nup62c92=importdata('Nup62_Nup58_Nup54_32_9_surf2.txt');
Nup62c92=Nup62c92.data;
% Import Nup62c run10
Nup62c101=importdata('Nup62_Nup58_Nup54_32_10_surf1.txt');
Nup62c101=Nup62c101.data;
Nup62c102=importdata('Nup62_Nup58_Nup54_32_10_surf2.txt');
Nup62c102=Nup62c102.data;
% Import Nup62c run11
Nup62c111=importdata('Nup62_Nup58_Nup54_32_11_surf1.txt');
Nup62c111=Nup62c111.data;
Nup62c112=importdata('Nup62_Nup58_Nup54_32_11_surf2.txt');
Nup62c112=Nup62c112.data;
% Import Nup62c run12
Nup62c121=importdata('Nup62_Nup58_Nup54_32_12_surf1.txt');
Nup62c121=Nup62c121.data;
Nup62c122=importdata('Nup62_Nup58_Nup54_32_12_surf2.txt');
Nup62c122=Nup62c122.data;
% Import Nup62c run13
Nup62c131=importdata('Nup62_Nup58_Nup54_32_13_surf1.txt');
Nup62c131=Nup62c131.data;
Nup62c132=importdata('Nup62_Nup58_Nup54_32_13_surf2.txt');
Nup62c132=Nup62c132.data;
% Import Nup62c run14
Nup62c141=importdata('Nup62_Nup58_Nup54_32_14_surf1.txt');
Nup62c141=Nup62c141.data;
Nup62c142=importdata('Nup62_Nup58_Nup54_32_14_surf2.txt');
Nup62c142=Nup62c142.data;
% Import Nup62c run15
Nup62c151=importdata('Nup62_Nup58_Nup54_32_15_surf1.txt');
Nup62c151=Nup62c151.data;
Nup62c152=importdata('Nup62_Nup58_Nup54_32_15_surf2.txt');
Nup62c152=Nup62c152.data;

t=[1:5000]/1000;
t2=[1:5001]/1000;

% Compare meassures on the same system
figure;
hold on;
plot(t,Nup62c11,'Linewidth',3,'Color',color1);
plot(t2,Nup62c21,'Linewidth',3,'Color',color2);
plot(t,Nup62c31,'Linewidth',3,'Color',color3);
plot(t,Nup62c41,'Linewidth',3,'Color',color4);
plot(t,Nup62c51,'Linewidth',3,'Color',color5);
plot(t,Nup62c61,'Linewidth',3,'Color',color6);
plot(t,Nup62c71,'Linewidth',3,'Color',color7);
plot(t,Nup62c81,'Linewidth',3,'Color',color8);
plot(t,Nup62c91,'--','Linewidth',3,'Color',color1);
plot(t,Nup62c101,'--','Linewidth',3,'Color',color2);
plot(t,Nup62c111,'--','Linewidth',3,'Color',color3);
plot(t,Nup62c121,'--','Linewidth',3,'Color',color4);
plot(t,Nup62c131,'--','Linewidth',3,'Color',color5);
plot(t,Nup62c141,'--','Linewidth',3,'Color',color6);
plot(t,Nup62c151,'--','Linewidth',3,'Color',color7);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4,'yscale','log');
% legend({'method 1','methdod 2'},'location','northwest','FontSize',52,'FontName','Helvetica');
axis([0 5 0 100]);
box on;

% arrays with simulations - remove outliers by radius (sim 4, 8, and 9) -----------------------------------------------------
surf_WT_tot1=[mean(Nup62c11);mean(Nup62c21);mean(Nup62c31);mean(Nup62c51);mean(Nup62c61);mean(Nup62c71);mean(Nup62c101);mean(Nup62c111);mean(Nup62c121);mean(Nup62c131);mean(Nup62c141);mean(Nup62c151)];
surf_WT_tot2=[mean(Nup62c12);mean(Nup62c22);mean(Nup62c32);mean(Nup62c52);mean(Nup62c62);mean(Nup62c72);mean(Nup62c102);mean(Nup62c112);mean(Nup62c122);mean(Nup62c132);mean(Nup62c142);mean(Nup62c152)];

mean(surf_WT_tot1)
mean(surf_WT_tot2)