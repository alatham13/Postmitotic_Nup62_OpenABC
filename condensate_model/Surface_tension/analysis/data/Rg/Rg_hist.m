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

% Constant volume sims with 32 copies -------------------------------------
Nup62c1=importdata('Nup62_Nup58_Nup54_32_Rg.txt',' ',1);
Nup62c1=Nup62c1.data;
Nup62c2=importdata('Nup62_Nup58_Nup54_32_2_Rg.txt',' ',1);
Nup62c2=Nup62c2.data;
Nup62c3=importdata('Nup62_Nup58_Nup54_32_3_Rg.txt',' ',1);
Nup62c3=Nup62c3.data;
Nup62c4=importdata('Nup62_Nup58_Nup54_32_4_Rg.txt',' ',1);
Nup62c4=Nup62c4.data;
Nup62c5=importdata('Nup62_Nup58_Nup54_32_5_Rg.txt',' ',1);
Nup62c5=Nup62c5.data;
Nup62c6=importdata('Nup62_Nup58_Nup54_32_6_Rg.txt',' ',1);
Nup62c6=Nup62c6.data;
Nup62c7=importdata('Nup62_Nup58_Nup54_32_7_Rg.txt',' ',1);
Nup62c7=Nup62c7.data;
Nup62c8=importdata('Nup62_Nup58_Nup54_32_8_Rg.txt',' ',1);
Nup62c8=Nup62c8.data;
Nup62c9=importdata('Nup62_Nup58_Nup54_32_9_Rg.txt',' ',1);
Nup62c9=Nup62c9.data;
Nup62c10=importdata('Nup62_Nup58_Nup54_32_10_Rg.txt',' ',1);
Nup62c10=Nup62c10.data;
Nup62c11=importdata('Nup62_Nup58_Nup54_32_11_Rg.txt',' ',1);
Nup62c11=Nup62c11.data;
Nup62c12=importdata('Nup62_Nup58_Nup54_32_12_Rg.txt',' ',1);
Nup62c12=Nup62c12.data;
Nup62c13=importdata('Nup62_Nup58_Nup54_32_13_Rg.txt',' ',1);
Nup62c13=Nup62c13.data;
Nup62c14=importdata('Nup62_Nup58_Nup54_32_14_Rg.txt',' ',1);
Nup62c14=Nup62c14.data;
Nup62c15=importdata('Nup62_Nup58_Nup54_32_15_Rg.txt',' ',1);
Nup62c15=Nup62c15.data;

% Monomer sims ------------------------------------------------------------
Nup62c_mono1=importdata('Nup62_Nup58_Nup54_monomer.txt',' ',1);
Nup62c_mono1=Nup62c_mono1.data;
Nup62c_mono2=importdata('Nup62_Nup58_Nup54_monomer2.txt',' ',1);
Nup62c_mono2=Nup62c_mono2.data;
Nup62c_mono3=importdata('Nup62_Nup58_Nup54_monomer3.txt',' ',1);
Nup62c_mono3=Nup62c_mono3.data;
Nup62c_mono4=importdata('Nup62_Nup58_Nup54_monomer4.txt',' ',1);
Nup62c_mono4=Nup62c_mono4.data;
Nup62c_mono5=importdata('Nup62_Nup58_Nup54_monomer5.txt',' ',1);
Nup62c_mono5=Nup62c_mono5.data;
Nup62c_mono=[Nup62c_mono1;Nup62c_mono2;Nup62c_mono3;Nup62c_mono4;Nup62c_mono5];

% Examine total. Remove simulations based on condensation
Nup62c_tot1=[Nup62c1;Nup62c3;Nup62c4;Nup62c5;Nup62c6;Nup62c7;Nup62c9;Nup62c10;Nup62c11;Nup62c12];

edges=[2.95:0.1:10.05];
[N1,edges]=histcounts(Nup62c_mono(:),edges,'Normalization','pdf');
[N2,edges]=histcounts(Nup62c_tot1(:),edges,'Normalization','pdf');
edges2=[3:0.1:10];

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(edges2,N1,'Linewidth',3,'Color',blue1);
plot(edges2,N2,'Linewidth',3,'Color',red1);
legend({' dilute',' condensate'},'location','northeast');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([3 10 0 1.5]);
legend boxoff;
box on;
print(fig,'Rg_comp.png','-dpng');

mean(Nup62c_mono)
mean(mean(Nup62c_tot1))

