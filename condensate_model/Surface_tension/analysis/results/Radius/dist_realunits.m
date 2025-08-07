% code to plot and 

clear all;
close all;

colorB=[0, 0.4470, 0.7410];
colorA=[0.8500, 0.3250, 0.0980];
color3=[0.9290, 0.6940, 0.1250];
color4=[0.4940, 0.1840, 0.5560];
color5=[0.4660, 0.6740, 0.1880];
color6=[0.3010, 0.7450, 0.9330];
color7=[0.6350, 0.0780, 0.1840];
color8=[160,  82,  45]/255;

blue1=[0,0,128]/255;
blue2=[65,105,225]/255;
blue3=[135,206,250]/255;

red1=[178,34,34]/255;
red2=[220,20,60]/255;
red3=[250,128,114]/255;

dZ=500/100;

% Nup62c -------------------------------------------------------------
cluster=importdata('Nup62_Nup58_Nup54_32_size.txt');
check_cluster=min(cluster(:,1))
Nup62c1=importdata('Nup62_Nup58_Nup54_32_fitted.txt',' ',1);
Nup62c1=Nup62c1.data;

cluster=importdata('Nup62_Nup58_Nup54_32_2_size.txt');
check_cluster=min(cluster(:,1))
Nup62c2=importdata('Nup62_Nup58_Nup54_32_2_fitted.txt',' ',1);
Nup62c2=Nup62c2.data;

cluster=importdata('Nup62_Nup58_Nup54_32_3_size.txt');
check_cluster=min(cluster(:,1))
Nup62c3=importdata('Nup62_Nup58_Nup54_32_3_fitted.txt',' ',1);
Nup62c3=Nup62c3.data;

cluster=importdata('Nup62_Nup58_Nup54_32_4_size.txt');
check_cluster=min(cluster(:,1))
Nup62c4=importdata('Nup62_Nup58_Nup54_32_4_fitted.txt',' ',1);
Nup62c4=Nup62c4.data;

cluster=importdata('Nup62_Nup58_Nup54_32_5_size.txt');
check_cluster=min(cluster(:,1))
Nup62c5=importdata('Nup62_Nup58_Nup54_32_5_fitted.txt',' ',1);
Nup62c5=Nup62c5.data;

cluster=importdata('Nup62_Nup58_Nup54_32_6_size.txt');
check_cluster=min(cluster(:,1))
Nup62c6=importdata('Nup62_Nup58_Nup54_32_6_fitted.txt',' ',1);
Nup62c6=Nup62c6.data;

cluster=importdata('Nup62_Nup58_Nup54_32_7_size.txt');
check_cluster=min(cluster(:,1))
Nup62c7=importdata('Nup62_Nup58_Nup54_32_7_fitted.txt',' ',1);
Nup62c7=Nup62c7.data;

cluster=importdata('Nup62_Nup58_Nup54_32_8_size.txt');
check_cluster=min(cluster(:,1))
Nup62c8=importdata('Nup62_Nup58_Nup54_32_8_fitted.txt',' ',1);
Nup62c8=Nup62c8.data;

cluster=importdata('Nup62_Nup58_Nup54_32_9_size.txt');
check_cluster=min(cluster(:,1))
Nup62c9=importdata('Nup62_Nup58_Nup54_32_9_fitted.txt',' ',1);
Nup62c9=Nup62c9.data;

cluster=importdata('Nup62_Nup58_Nup54_32_10_size.txt');
check_cluster=min(cluster(:,1))
Nup62c10=importdata('Nup62_Nup58_Nup54_32_10_fitted.txt',' ',1);
Nup62c10=Nup62c10.data;

cluster=importdata('Nup62_Nup58_Nup54_32_11_size.txt');
check_cluster=min(cluster(:,1))
Nup62c11=importdata('Nup62_Nup58_Nup54_32_11_fitted.txt',' ',1);
Nup62c11=Nup62c11.data;

cluster=importdata('Nup62_Nup58_Nup54_32_12_size.txt');
check_cluster=min(cluster(:,1))
Nup62c12=importdata('Nup62_Nup58_Nup54_32_12_fitted.txt',' ',1);
Nup62c12=Nup62c12.data;

cluster=importdata('Nup62_Nup58_Nup54_32_13_size.txt');
check_cluster=min(cluster(:,1))
Nup62c13=importdata('Nup62_Nup58_Nup54_32_13_fitted.txt',' ',1);
Nup62c13=Nup62c13.data;

cluster=importdata('Nup62_Nup58_Nup54_32_14_size.txt');
check_cluster=min(cluster(:,1))
Nup62c14=importdata('Nup62_Nup58_Nup54_32_14_fitted.txt',' ',1);
Nup62c14=Nup62c14.data;

cluster=importdata('Nup62_Nup58_Nup54_32_15_size.txt');
check_cluster=min(cluster(:,1))
Nup62c15=importdata('Nup62_Nup58_Nup54_32_15_fitted.txt',' ',1);
Nup62c15=Nup62c15.data;

% Calculate radius and check for outliers. Removing samples more than 3
% median absolute deviations from the median
R=[14.180254977001598;15.359661208680578;14.853634802839506;17.066738266284624;7.682436834512304;17.501196971375958;15.731348592877895;3.818679466059395;999.8870668213114;15.754040879893715;13.735494835346959;14.077873616776097;10.341160836276845;17.006827799398955;16.954927542663366];

out=isoutlier(R)

R2=R(out~=1)
mean(R2)

% Remove simulation: 4, 8, 9

Nup62c_mean=mean([Nup62c1(:,2),Nup62c2(:,2),Nup62c3(:,2),Nup62c5(:,2),Nup62c6(:,2),Nup62c7(:,2),Nup62c10(:,2),Nup62c11(:,2),Nup62c12(:,2),Nup62c13(:,2),Nup62c14(:,2),Nup62c15(:,2)],2);
Nup62c_std=std([Nup62c1(:,2),Nup62c2(:,2),Nup62c3(:,2),Nup62c5(:,2),Nup62c6(:,2),Nup62c7(:,2),Nup62c10(:,2),Nup62c11(:,2),Nup62c12(:,2),Nup62c13(:,2),Nup62c14(:,2),Nup62c15(:,2)],0,2);

Nup62c_fit_mean=mean([Nup62c1(:,3),Nup62c2(:,3),Nup62c3(:,3),Nup62c5(:,3),Nup62c6(:,3),Nup62c7(:,3),Nup62c10(:,3),Nup62c11(:,3),Nup62c12(:,3),Nup62c13(:,3),Nup62c14(:,3),Nup62c15(:,3)],2);
Nup62c_fit_std=std([Nup62c1(:,3),Nup62c2(:,3),Nup62c3(:,3),Nup62c5(:,3),Nup62c6(:,3),Nup62c7(:,3),Nup62c10(:,3),Nup62c11(:,3),Nup62c12(:,3),Nup62c13(:,3),Nup62c14(:,3),Nup62c15(:,3)],0,2);


figure;
hold on;
plot(Nup62c1(:,1),Nup62c_mean,'Linewidth',6,'Color',colorB);
plot(Nup62c1(:,1),Nup62c_fit_mean,'Linewidth',6,'Color',colorA);
errorbar(Nup62c1(:,1),Nup62c_mean,Nup62c_std,'Linewidth',6,'Color',colorB);
errorbar(Nup62c1(:,1),Nup62c_fit_mean,Nup62c_fit_std,'Linewidth',6,'Color',colorA);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4);
legend({'Simulation','Fitted'},'location','northeast','FontSize',52,'FontName','Helvetica');
axis([0 max(Nup62c1(:,1)) 0 600]);
box on;

X=[mean(R2);mean(R2)];
Y=[0;600];

figure;
hold on;
plot(Nup62c1(:,1),Nup62c_mean,'Linewidth',6,'Color',colorB);
plot(X,Y,'k--','Linewidth',6);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4);
% legend({'Simulation','Fitted'},'location','northeast','FontSize',52,'FontName','Helvetica');
axis([0 max(Nup62c1(:,1)) 0 600]);
box on;

V=(4/3)*pi*mean(R2)^3

