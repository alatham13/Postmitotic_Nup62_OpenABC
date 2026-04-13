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
color9=[0.6350, 0.0780, 0.1840];

color1=[0, 0.4470, 0.7410];
color2=[0.8500, 0.3250, 0.0980];
color3=[0.9290, 0.6940, 0.1250];
color4=[0.4940, 0.1840, 0.5560];
color5=[0.4660, 0.6740, 0.1880];
color6=[0.3010, 0.7450, 0.9330];

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
% Convert mg / mL to M
mass=169537.1900000013; % mass in g / mol
Nup62c(:,2)=Nup62c(:,2)./mass;
% convert M to micro M
Nup62c(:,2)=Nup62c(:,2).*10^6;

% Nup62c09 -------------------------------------------------------------
Nup62c09=importdata('Nup62_Nup58_Nup54_eps0.9_hist.txt',' ',1);
Nup62c09=Nup62c09.data;
% Normalize
Nup62c09(:,2)=Nup62c09(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
xy=36.22019769301058;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
Nup62c09(:,2)=Nup62c09(:,2)*1.6605*10^(-21);
Nup62c09(:,2)=Nup62c09(:,2)./(dV);
% Convert mg / mL to M
mass=169537.1900000013; % mass in g / mol
Nup62c09(:,2)=Nup62c09(:,2)./mass;
% convert M to micro M
Nup62c09(:,2)=Nup62c09(:,2).*10^6;

% Nup62c07 -------------------------------------------------------------
Nup62c07=importdata('Nup62_Nup58_Nup54_eps0.7_hist.txt',' ',1);
Nup62c07=Nup62c07.data;
% Normalize
Nup62c07(:,2)=Nup62c07(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
xy=38.04779485532567;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
Nup62c07(:,2)=Nup62c07(:,2)*1.6605*10^(-21);
Nup62c07(:,2)=Nup62c07(:,2)./(dV);
% Convert mg / mL to M
mass=169537.1900000013; % mass in g / mol
Nup62c07(:,2)=Nup62c07(:,2)./mass;
% convert M to micro M
Nup62c07(:,2)=Nup62c07(:,2).*10^6;

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(Nup62c(:,1)/10,Nup62c(:,2),'Linewidth',3,'Color',color1);
plot(Nup62c09(:,1)/10,Nup62c09(:,2),'Linewidth',3,'Color',color2);
plot(Nup62c07(:,1)/10,Nup62c07(:,2),'Linewidth',3,'Color',color5);
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
legend({' 0%',' 1.5%',' 5.0%'},'location','northeast','FontSize',36,'FontName','Helvetica');
axis([-250 250 0 3000]);
legend boxoff;
box on;
print(fig,'eps_mM.png','-dpng');

exp=[2;2];
edges=[-250;250];

y=[10^(-4);10^5;10^5;10^(-4)];
x1=[Nup62c09(100);Nup62c09(100);Nup62c09(91);Nup62c09(91)]/10;
x2=[Nup62c09(1);Nup62c09(1);Nup62c09(10);Nup62c09(10)]/10;

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(Nup62c(:,1)/10,Nup62c(:,2),'Linewidth',3,'Color',color1);
plot(Nup62c09(:,1)/10,Nup62c09(:,2),'Linewidth',3,'Color',color2);
plot(Nup62c07(:,1)/10,Nup62c07(:,2),'Linewidth',3,'Color',color5);
plot(edges,exp,'k--','Linewidth',3);
fill(x1,y,'y', 'FaceAlpha', 0.3,'LineStyle','none');
fill(x2,y,'y', 'FaceAlpha', 0.3,'LineStyle','none');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
legend({' 0%',' 1.5%',' 5.0%'},'location','north','FontSize',36,'FontName','Helvetica','NumColumns',2);
yscale log;
axis([-250 250 10^(-1) 10^5]);
legend boxoff;
box on;
print(fig,'eps_log_mM.png','-dpng');

csat=mean([Nup62c09(91:100,2);Nup62c09(1:10,2)])

