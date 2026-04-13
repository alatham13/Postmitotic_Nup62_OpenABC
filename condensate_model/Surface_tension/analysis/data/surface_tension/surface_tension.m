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

function ci = conf_int(dist)
    xbar = mean(dist);
    n = max(size(dist));
    se = std(dist)/sqrt(n);
    crit=[1.96,-1.96];
    ci = xbar + crit*se;
end


% Calculate Pave ---------------------------------------------------------
Pave1=importdata('Nup62_Nup58_Nup54_32_surf_Pave_v2.txt');
Pave1=Pave1.data;
Pave2=importdata('Nup62_Nup58_Nup54_32_2_surf_Pave_v2.txt');
Pave2=Pave2.data;
Pave3=importdata('Nup62_Nup58_Nup54_32_3_surf_Pave_v2.txt');
Pave3=Pave3.data;
Pave4=importdata('Nup62_Nup58_Nup54_32_4_surf_Pave_v2.txt');
Pave4=Pave4.data;
Pave5=importdata('Nup62_Nup58_Nup54_32_5_surf_Pave_v2.txt');
Pave5=Pave5.data;
Pave6=importdata('Nup62_Nup58_Nup54_32_6_surf_Pave_v2.txt');
Pave6=Pave6.data;
Pave7=importdata('Nup62_Nup58_Nup54_32_7_surf_Pave_v2.txt');
Pave7=Pave7.data;
Pave8=importdata('Nup62_Nup58_Nup54_32_8_surf_Pave_v2.txt');
Pave8=Pave8.data;
Pave9=importdata('Nup62_Nup58_Nup54_32_9_surf_Pave_v2.txt');
Pave9=Pave9.data;
Pave10=importdata('Nup62_Nup58_Nup54_32_10_surf_Pave_v2.txt');
Pave10=Pave10.data;
Pave11=importdata('Nup62_Nup58_Nup54_32_11_surf_Pave_v2.txt');
Pave11=Pave11.data;
Pave12=importdata('Nup62_Nup58_Nup54_32_12_surf_Pave_v2.txt');
Pave12=Pave12.data;
Pave13=importdata('Nup62_Nup58_Nup54_32_13_surf_Pave_v2.txt');
Pave13=Pave13.data;
Pave14=importdata('Nup62_Nup58_Nup54_32_14_surf_Pave_v2.txt');
Pave14=Pave14.data;
Pave15=importdata('Nup62_Nup58_Nup54_32_15_surf_Pave_v2.txt');
Pave15=Pave15.data;

% Remove only simulations with multiple condensates
surf_WT_Pave=[mean(Pave1);mean(Pave3);mean(Pave4);mean(Pave5);mean(Pave6);mean(Pave7);mean(Pave9);mean(Pave10);mean(Pave11);mean(Pave12)];
mean(surf_WT_Pave)
ci_Pave=conf_int(surf_WT_Pave)
