clear all;
close all;
step=3000; % number of timesteps

% color1_v1=[255, 0, 255]/255; % Nup214_v1
% color1_v2=[236, 0, 140]/255; % Nup214_v2
% color1=[241, 64, 169]/255; % Nup62c
% color2=[245, 128, 198]/255; % Nup214c
% color3=[188, 166, 233]/255; % Nup54
% color4=[228, 219, 246]/255; % Nup58_CT
% color5=[110, 84, 164]/255; % Nup58_NT
% color6=[147, 112, 219]/255; % Nup62

color7=[102, 45, 145]/255; % Nup62c
color8_v1=[255, 0, 255]/255; % Nup214c_v1
color8_v2=[236, 0, 140]/255; % Nup214c_v2

dZ=500/100;

% Nup62c -------------------------------------------------------------
Nup62c=importdata('Nup62_Nup58_Nup54_hist.txt',' ',1);
Nup62c=Nup62c.data;
% Normalize
Nup62c(:,2)=Nup62c(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
xy=35.7133310893075;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
Nup62c(:,2)=Nup62c(:,2)*1.6605*10^(-21);
Nup62c(:,2)=Nup62c(:,2)./(dV);

figure;
hold on;
plot(Nup62c(:,1)/10,Nup62c(:,2),'Linewidth',6,'Color',color7);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4);
legend({'Nup62-Nup58-Nup54'},'location','northeast','FontSize',52,'FontName','Helvetica');
axis([-250 250 0 650]);
box on;

csat=zeros(1,1);
csat(1)=mean([Nup62c(1:10,2);Nup62c(91:100,2)]);

cdense=zeros(1,1);
cdense(1)=mean(Nup62c(49:52,2));
