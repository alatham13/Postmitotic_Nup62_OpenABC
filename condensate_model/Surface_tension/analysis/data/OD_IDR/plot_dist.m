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

Nup62c_color=[102, 45, 145]/255; % Nup62c
grey=[0.35, 0.35, 0.35];

blue1=[0,0,128]/255;
blue2=[65,105,225]/255;
blue3=[135,206,250]/255;

red1=[178,34,34]/255;
red2=[220,20,60]/255;
red3=[250,128,114]/255;

dZ=500/100;

function normalized_density = NormalizeDensity(density_data,cluster_size_fn)
    radius=density_data(:,1)/10;
    density=density_data(:,2);
    % % Correct radius to start / end at the middle of the bin instead of the
    % % beginning / end. (was adjusted for plotting)
    dR = radius(3) - radius(2);
    % radius(1)=radius(1)+dR/2;
    % radius(max(size(radius)))=radius(max(size(radius)))-dR/2;
    % Normalize density
    cluster_size=importdata(cluster_size_fn);
    sim_length=max(size(cluster_size));
    normalized_density=density/sim_length;
    for i = 1:max(size(normalized_density))
        index0=i-1;
        r1=(index0+1)*dR;
        V1=(4/3) * pi * (r1 ^ 3);
        r2 = (index0) * dR;
        V2 = (4 / 3) * pi * (r2 ^ 3);
        dV = (V1-V2)/ (10^(21));
        normalized_density(i)=normalized_density(i)*1.6605*10^(-21);
        normalized_density(i)=normalized_density(i)/dV;
    end
end

% Import and normalize data from simulation 1 -----------------------------
Nup62c_clusterfn1='Nup62_Nup58_Nup54_32_size.txt';
Nup62c1=importdata('Nup62_Nup58_Nup54_32_total.txt',' ',1);
Nup62c1=Nup62c1.data;
Nup62c1_norm=NormalizeDensity(Nup62c1,Nup62c_clusterfn1);
OD1=importdata('Nup62_Nup58_Nup54_32_OD.txt',' ',1);
OD1=OD1.data;
OD1_norm=NormalizeDensity(OD1,Nup62c_clusterfn1);
IDR1=importdata('Nup62_Nup58_Nup54_32_IDR.txt',' ',1);
IDR1=IDR1.data;
IDR1_norm=NormalizeDensity(IDR1,Nup62c_clusterfn1);

% Import and normalize data from simulation 2 -----------------------------
Nup62c_clusterfn2='Nup62_Nup58_Nup54_32_2_size.txt';
Nup62c2=importdata('Nup62_Nup58_Nup54_32_2_total.txt',' ',1);
Nup62c2=Nup62c2.data;
Nup62c2_norm=NormalizeDensity(Nup62c2,Nup62c_clusterfn2);
OD2=importdata('Nup62_Nup58_Nup54_32_2_OD.txt',' ',1);
OD2=OD2.data;
OD2_norm=NormalizeDensity(OD2,Nup62c_clusterfn2);
IDR2=importdata('Nup62_Nup58_Nup54_32_2_IDR.txt',' ',1);
IDR2=IDR2.data;
IDR2_norm=NormalizeDensity(IDR2,Nup62c_clusterfn2);

% Import and normalize data from simulation 3 -----------------------------
Nup62c_clusterfn3='Nup62_Nup58_Nup54_32_3_size.txt';
Nup62c3=importdata('Nup62_Nup58_Nup54_32_3_total.txt',' ',1);
Nup62c3=Nup62c3.data;
Nup62c3_norm=NormalizeDensity(Nup62c3,Nup62c_clusterfn3);
OD3=importdata('Nup62_Nup58_Nup54_32_3_OD.txt',' ',1);
OD3=OD3.data;
OD3_norm=NormalizeDensity(OD3,Nup62c_clusterfn3);
IDR3=importdata('Nup62_Nup58_Nup54_32_3_IDR.txt',' ',1);
IDR3=IDR3.data;
IDR3_norm=NormalizeDensity(IDR3,Nup62c_clusterfn3);

% Import and normalize data from simulation 4 -----------------------------
Nup62c_clusterfn4='Nup62_Nup58_Nup54_32_4_size.txt';
Nup62c4=importdata('Nup62_Nup58_Nup54_32_4_total.txt',' ',1);
Nup62c4=Nup62c4.data;
Nup62c4_norm=NormalizeDensity(Nup62c4,Nup62c_clusterfn4);
OD4=importdata('Nup62_Nup58_Nup54_32_4_OD.txt',' ',1);
OD4=OD4.data;
OD4_norm=NormalizeDensity(OD4,Nup62c_clusterfn4);
IDR4=importdata('Nup62_Nup58_Nup54_32_4_IDR.txt',' ',1);
IDR4=IDR4.data;
IDR4_norm=NormalizeDensity(IDR4,Nup62c_clusterfn4);

% Import and normalize data from simulation 5 -----------------------------
Nup62c_clusterfn5='Nup62_Nup58_Nup54_32_5_size.txt';
Nup62c5=importdata('Nup62_Nup58_Nup54_32_5_total.txt',' ',1);
Nup62c5=Nup62c5.data;
Nup62c5_norm=NormalizeDensity(Nup62c5,Nup62c_clusterfn5);
OD5=importdata('Nup62_Nup58_Nup54_32_5_OD.txt',' ',1);
OD5=OD5.data;
OD5_norm=NormalizeDensity(OD5,Nup62c_clusterfn5);
IDR5=importdata('Nup62_Nup58_Nup54_32_5_IDR.txt',' ',1);
IDR5=IDR5.data;
IDR5_norm=NormalizeDensity(IDR5,Nup62c_clusterfn5);

% Import and normalize data from simulation 6 -----------------------------
Nup62c_clusterfn6='Nup62_Nup58_Nup54_32_6_size.txt';
Nup62c6=importdata('Nup62_Nup58_Nup54_32_6_total.txt',' ',1);
Nup62c6=Nup62c6.data;
Nup62c6_norm=NormalizeDensity(Nup62c6,Nup62c_clusterfn6);
OD6=importdata('Nup62_Nup58_Nup54_32_6_OD.txt',' ',1);
OD6=OD6.data;
OD6_norm=NormalizeDensity(OD6,Nup62c_clusterfn6);
IDR6=importdata('Nup62_Nup58_Nup54_32_6_IDR.txt',' ',1);
IDR6=IDR6.data;
IDR6_norm=NormalizeDensity(IDR6,Nup62c_clusterfn6);

% Import and normalize data from simulation 7 -----------------------------
Nup62c_clusterfn7='Nup62_Nup58_Nup54_32_7_size.txt';
Nup62c7=importdata('Nup62_Nup58_Nup54_32_7_total.txt',' ',1);
Nup62c7=Nup62c7.data;
Nup62c7_norm=NormalizeDensity(Nup62c7,Nup62c_clusterfn7);
OD7=importdata('Nup62_Nup58_Nup54_32_7_OD.txt',' ',1);
OD7=OD7.data;
OD7_norm=NormalizeDensity(OD7,Nup62c_clusterfn7);
IDR7=importdata('Nup62_Nup58_Nup54_32_7_IDR.txt',' ',1);
IDR7=IDR7.data;
IDR7_norm=NormalizeDensity(IDR7,Nup62c_clusterfn7);

% Import and normalize data from simulation 8 -----------------------------
Nup62c_clusterfn8='Nup62_Nup58_Nup54_32_8_size.txt';
Nup62c8=importdata('Nup62_Nup58_Nup54_32_8_total.txt',' ',1);
Nup62c8=Nup62c8.data;
Nup62c8_norm=NormalizeDensity(Nup62c8,Nup62c_clusterfn8);
OD8=importdata('Nup62_Nup58_Nup54_32_8_OD.txt',' ',1);
OD8=OD8.data;
OD8_norm=NormalizeDensity(OD8,Nup62c_clusterfn8);
IDR8=importdata('Nup62_Nup58_Nup54_32_8_IDR.txt',' ',1);
IDR8=IDR8.data;
IDR8_norm=NormalizeDensity(IDR8,Nup62c_clusterfn8);

% Import and normalize data from simulation 9 -----------------------------
Nup62c_clusterfn9='Nup62_Nup58_Nup54_32_9_size.txt';
Nup62c9=importdata('Nup62_Nup58_Nup54_32_9_total.txt',' ',1);
Nup62c9=Nup62c9.data;
Nup62c9_norm=NormalizeDensity(Nup62c9,Nup62c_clusterfn9);
OD9=importdata('Nup62_Nup58_Nup54_32_9_OD.txt',' ',1);
OD9=OD9.data;
OD9_norm=NormalizeDensity(OD9,Nup62c_clusterfn9);
IDR9=importdata('Nup62_Nup58_Nup54_32_9_IDR.txt',' ',1);
IDR9=IDR9.data;
IDR9_norm=NormalizeDensity(IDR9,Nup62c_clusterfn9);

% Import and normalize data from simulation 10 ----------------------------
Nup62c_clusterfn10='Nup62_Nup58_Nup54_32_10_size.txt';
Nup62c10=importdata('Nup62_Nup58_Nup54_32_10_total.txt',' ',1);
Nup62c10=Nup62c10.data;
Nup62c10_norm=NormalizeDensity(Nup62c10,Nup62c_clusterfn10);
OD10=importdata('Nup62_Nup58_Nup54_32_10_OD.txt',' ',1);
OD10=OD10.data;
OD10_norm=NormalizeDensity(OD10,Nup62c_clusterfn10);
IDR10=importdata('Nup62_Nup58_Nup54_32_10_IDR.txt',' ',1);
IDR10=IDR10.data;
IDR10_norm=NormalizeDensity(IDR10,Nup62c_clusterfn10);

% Import and normalize data from simulation 11 ----------------------------
Nup62c_clusterfn11='Nup62_Nup58_Nup54_32_11_size.txt';
Nup62c11=importdata('Nup62_Nup58_Nup54_32_11_total.txt',' ',1);
Nup62c11=Nup62c11.data;
Nup62c11_norm=NormalizeDensity(Nup62c11,Nup62c_clusterfn11);
OD11=importdata('Nup62_Nup58_Nup54_32_11_OD.txt',' ',1);
OD11=OD11.data;
OD11_norm=NormalizeDensity(OD11,Nup62c_clusterfn11);
IDR11=importdata('Nup62_Nup58_Nup54_32_11_IDR.txt',' ',1);
IDR11=IDR11.data;
IDR11_norm=NormalizeDensity(IDR11,Nup62c_clusterfn11);

% Import and normalize data from simulation 12 ----------------------------
Nup62c_clusterfn12='Nup62_Nup58_Nup54_32_12_size.txt';
Nup62c12=importdata('Nup62_Nup58_Nup54_32_12_total.txt',' ',1);
Nup62c12=Nup62c12.data;
Nup62c12_norm=NormalizeDensity(Nup62c12,Nup62c_clusterfn12);
OD12=importdata('Nup62_Nup58_Nup54_32_12_OD.txt',' ',1);
OD12=OD12.data;
OD12_norm=NormalizeDensity(OD12,Nup62c_clusterfn12);
IDR12=importdata('Nup62_Nup58_Nup54_32_12_IDR.txt',' ',1);
IDR12=IDR12.data;
IDR12_norm=NormalizeDensity(IDR12,Nup62c_clusterfn12);

% Import and normalize data from simulation 13 ----------------------------
Nup62c_clusterfn13='Nup62_Nup58_Nup54_32_13_size.txt';
Nup62c13=importdata('Nup62_Nup58_Nup54_32_13_total.txt',' ',1);
Nup62c13=Nup62c13.data;
Nup62c13_norm=NormalizeDensity(Nup62c13,Nup62c_clusterfn13);
OD13=importdata('Nup62_Nup58_Nup54_32_13_OD.txt',' ',1);
OD13=OD13.data;
OD13_norm=NormalizeDensity(OD13,Nup62c_clusterfn13);
IDR13=importdata('Nup62_Nup58_Nup54_32_13_IDR.txt',' ',1);
IDR13=IDR13.data;
IDR13_norm=NormalizeDensity(IDR13,Nup62c_clusterfn13);

% Import and normalize data from simulation 14 ----------------------------
Nup62c_clusterfn14='Nup62_Nup58_Nup54_32_14_size.txt';
Nup62c14=importdata('Nup62_Nup58_Nup54_32_14_total.txt',' ',1);
Nup62c14=Nup62c14.data;
Nup62c14_norm=NormalizeDensity(Nup62c14,Nup62c_clusterfn14);
OD14=importdata('Nup62_Nup58_Nup54_32_14_OD.txt',' ',1);
OD14=OD14.data;
OD14_norm=NormalizeDensity(OD14,Nup62c_clusterfn14);
IDR14=importdata('Nup62_Nup58_Nup54_32_14_IDR.txt',' ',1);
IDR14=IDR14.data;
IDR14_norm=NormalizeDensity(IDR14,Nup62c_clusterfn14);

% Import and normalize data from simulation 15 ----------------------------
Nup62c_clusterfn15='Nup62_Nup58_Nup54_32_15_size.txt';
Nup62c15=importdata('Nup62_Nup58_Nup54_32_15_total.txt',' ',1);
Nup62c15=Nup62c15.data;
Nup62c15_norm=NormalizeDensity(Nup62c15,Nup62c_clusterfn15);
OD15=importdata('Nup62_Nup58_Nup54_32_15_OD.txt',' ',1);
OD15=OD15.data;
OD15_norm=NormalizeDensity(OD15,Nup62c_clusterfn15);
IDR15=importdata('Nup62_Nup58_Nup54_32_15_IDR.txt',' ',1);
IDR15=IDR15.data;
IDR15_norm=NormalizeDensity(IDR15,Nup62c_clusterfn15);

% calculate with 1,3,4,5,6,7,9,10,11,12 --------------------------------------
OD_centered_tot=[Nup62c1_norm,Nup62c4_norm,Nup62c6_norm,Nup62c9_norm,Nup62c11_norm,Nup62c12_norm];
OD_centered_OD=[OD1_norm,OD4_norm,OD6_norm,OD9_norm,OD11_norm,OD12_norm];
OD_centered_IDR=[IDR1_norm,IDR4_norm,IDR6_norm,IDR9_norm,IDR11_norm,IDR12_norm];

% calculate with 1,3,4,5,6,7,9,10,11,12 --------------------------------------
IDR_centered_tot=[Nup62c3_norm,Nup62c5_norm,Nup62c7_norm,Nup62c10_norm];
IDR_centered_OD=[OD3_norm,OD5_norm,OD7_norm,OD10_norm];
IDR_centered_IDR=[IDR3_norm,IDR5_norm,IDR7_norm,IDR10_norm];

% plote results --------------------------------------------------------------

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(Nup62c1(:,1)./10,mean(OD_centered_tot,2),'Linewidth',3,'Color',color1);
plot(Nup62c1(:,1)./10,mean(OD_centered_OD,2),'Linewidth',3,'Color',Nup62c_color);
plot(Nup62c1(:,1)./10,mean(OD_centered_IDR,2),'Linewidth',3,'Color',grey);
legend({' total',' OD',' IDR'},'location','northeast');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([0 40 0 550]);
legend boxoff;
box on;
print(fig,'OD_centered_ave.png','-dpng');

% fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
% hold on;
% plot(Nup62c1(:,1)./10,OD_centered_tot,'Linewidth',3,'Color',color1);
% plot(Nup62c1(:,1)./10,OD_centered_OD,'Linewidth',3,'Color',Nup62c_color);
% plot(Nup62c1(:,1)./10,OD_centered_IDR,'Linewidth',3,'Color',grey);
% legend({' total',' OD',' IDR'},'location','northeast');
% set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
% axis([0 40 0 550]);
% legend boxoff;
% box on;
% print(fig,'OD_centered_tot.png','-dpng');

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(Nup62c1(:,1)./10,mean(IDR_centered_tot,2),'Linewidth',3,'Color',color1);
plot(Nup62c1(:,1)./10,mean(IDR_centered_OD,2),'Linewidth',3,'Color',Nup62c_color);
plot(Nup62c1(:,1)./10,mean(IDR_centered_IDR,2),'Linewidth',3,'Color',grey);
legend({' total',' OD',' IDR'},'location','northeast');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([0 40 0 550]);
legend boxoff;
box on;
print(fig,'IDR_centered_ave.png','-dpng');

% fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
% hold on;
% plot(Nup62c1(:,1)./10,IDR_centered_tot,'Linewidth',3,'Color',color1);
% plot(Nup62c1(:,1)./10,IDR_centered_OD,'Linewidth',3,'Color',Nup62c_color);
% plot(Nup62c1(:,1)./10,IDR_centered_IDR,'Linewidth',3,'Color',grey);
% legend({' total',' OD',' IDR'},'location','northeast');
% set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
% axis([0 40 0 550]);
% legend boxoff;
% box on;
% print(fig,'IDR_centered_tot.png','-dpng');

