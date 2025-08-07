clear all;
close all;
step=3000; % number of timesteps

% color1_v1=[255, 0, 255]/255; % Nup214_v1
% color1_v2=[236, 0, 140]/255; % Nup214_v2
color1=[241, 64, 169]/255; % Nup214_CT
color2=[245, 128, 198]/255; % Nup214_NT
color3=[188, 166, 233]/255; % Nup54
color4=[228, 219, 246]/255; % Nup58_CT
color5=[110, 84, 164]/255; % Nup58_NT
color6=[147, 112, 219]/255; % Nup62

% color7=[102, 45, 145]/255; % Nup62c
% color8_v1=[255, 0, 255]/255; % Nup214c_v1
% color8_v2=[236, 0, 140]/255; % Nup214c_v2

dZ=500/100;

% Nup214_CT -------------------------------------------------------------
Nup214_CT=importdata('Nup214_CT_hist.txt',' ',1);
Nup214_CT=Nup214_CT.data;

% Normalize
Nup214_CT(:,2)=Nup214_CT(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
xy=31.71072346748188;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
Nup214_CT(:,2)=Nup214_CT(:,2)*1.6605*10^(-21);
Nup214_CT(:,2)=Nup214_CT(:,2)./(dV);

% Nup214_NT -------------------------------------------------------------
Nup214_NT=importdata('Nup214_NT_hist.txt',' ',1);
Nup214_NT=Nup214_NT.data;

% Normalize
Nup214_NT(:,2)=Nup214_NT(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
xy=30.01076329927889;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
Nup214_NT(:,2)=Nup214_NT(:,2)*1.6605*10^(-21);
Nup214_NT(:,2)=Nup214_NT(:,2)./(dV);

% Nup54 -------------------------------------------------------------
Nup54=importdata('Nup54_hist.txt',' ',1);
Nup54=Nup54.data;

% Normalize
Nup54(:,2)=Nup54(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
xy=12.251833057173107;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
Nup54(:,2)=Nup54(:,2)*1.6605*10^(-21);
Nup54(:,2)=Nup54(:,2)./(dV);

% Nup58_CT -------------------------------------------------------------
Nup58_CT=importdata('Nup58_CT_hist.txt',' ',1);
Nup58_CT=Nup58_CT.data;

% Normalize
Nup58_CT(:,2)=Nup58_CT(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
xy=16.18014181486587;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
Nup58_CT(:,2)=Nup58_CT(:,2)*1.6605*10^(-21);
Nup58_CT(:,2)=Nup58_CT(:,2)./(dV);

% Nup58_NT -------------------------------------------------------------
Nup58_NT=importdata('Nup58_NT_hist.txt',' ',1);
Nup58_NT=Nup58_NT.data;

% Normalize
Nup58_NT(:,2)=Nup58_NT(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
xy=17.51341036795967;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
Nup58_NT(:,2)=Nup58_NT(:,2)*1.6605*10^(-21);
Nup58_NT(:,2)=Nup58_NT(:,2)./(dV);

% Nup62 -------------------------------------------------------------
Nup62=importdata('Nup62_hist.txt',' ',1);
Nup62=Nup62.data;

% Normalize
Nup62(:,2)=Nup62(:,2)./(step);
% Normalize the mass per step by the volume of each bin (in mL)
xy=19.31083753068954;
dV=dZ*xy*xy/(10^(21));
% Convert dalton to mg
Nup62(:,2)=Nup62(:,2)*1.6605*10^(-21);
Nup62(:,2)=Nup62(:,2)./(dV);

figure;
hold on;
plot(Nup214_CT(:,1)/10,Nup214_CT(:,2),'Linewidth',6,'Color',color1);
plot(Nup214_NT(:,1)/10,Nup214_NT(:,2),'Linewidth',6,'Color',color2);
plot(Nup54(:,1)/10,Nup54(:,2),'Linewidth',6,'Color',color3);
plot(Nup58_CT(:,1)/10,Nup58_CT(:,2),'Linewidth',6,'Color',color4);
plot(Nup58_NT(:,1)/10,Nup58_NT(:,2),'Linewidth',6,'Color',color5);
plot(Nup62(:,1)/10,Nup62(:,2),'Linewidth',6,'Color',color6);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4);
legend({'Nup214 CT','Nup214 NT','Nup54','Nup58 CT','Nup58 NT','Nup62'},'location','northeast','FontSize',52,'FontName','Helvetica');
axis([-250 250 0 650]);
box on;

csat=zeros(6,1);
csat(1)=mean([Nup214_CT(1:10,2);Nup214_CT(91:100,2)]);
csat(2)=mean([Nup214_NT(1:10,2);Nup214_NT(91:100,2)]);
csat(3)=mean([Nup54(1:10,2);Nup54(91:100,2)]);
csat(4)=mean([Nup58_CT(1:10,2);Nup58_CT(91:100,2)]);
csat(5)=mean([Nup58_NT(1:10,2);Nup58_NT(91:100,2)]);
csat(6)=mean([Nup62(1:10,2);Nup62(91:100,2)]);

cdense=zeros(6,1);
cdense(1)=mean(Nup214_CT(49:52,2));
cdense(2)=mean(Nup214_NT(49:52,2));
cdense(3)=mean(Nup54(49:52,2));
cdense(4)=mean(Nup58_CT(49:52,2));
cdense(5)=mean(Nup58_NT(49:52,2));
cdense(6)=mean(Nup62(49:52,2));
comp_c=cdense./csat;