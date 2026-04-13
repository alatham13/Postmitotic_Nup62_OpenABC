% code to plot and 

clear all;
close all;

color1=[0, 0.4470, 0.7410];
color2=[0.8500, 0.3250, 0.0980];
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

% Import Rg data for N49 --------------------------------------------------
N49_HDX0=importdata('N49_HDX-0_Rg.txt',' ',1);
N49_HDX0=N49_HDX0.data;
N49_HDX15=importdata('N49_HDX-1.5_Rg.txt',' ',1);
N49_HDX15=N49_HDX15.data;
N49_HDX5=importdata('N49_HDX-5.0_Rg.txt',' ',1);
N49_HDX5=N49_HDX5.data;

edges=[0.45:0.1:10.05];
[N1,edges]=histcounts(N49_HDX0(:),edges,'Normalization','pdf');
[N2,edges]=histcounts(N49_HDX15(:),edges,'Normalization','pdf');
[N3,edges]=histcounts(N49_HDX5(:),edges,'Normalization','pdf');
edges2=[0.5:0.1:10];

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(edges2,N1,'Linewidth',3,'Color',color5);
plot(edges2,N2,'Linewidth',3,'Color',color6);
plot(edges2,N3,'Linewidth',3,'Color',color7);
legend({' 0% HDX',' 1.5% HDX',' 5% HDX'},'location','northeast');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([0.5 3 0 2.0]);
legend boxoff;
box on;
print(fig,'N49_Rg.png','-dpng');


shift1=mean(N49_HDX15)-mean(N49_HDX0)

shift2=mean(N49_HDX5)-mean(N49_HDX0)
