
clear all;
close all;

color1=[0, 0.4470, 0.7410];
color2=[0.8500, 0.3250, 0.0980];
color3=[0.9290, 0.6940, 0.1250];
color4=[0.4940, 0.1840, 0.5560];
color5=[0.4660, 0.6740, 0.1880];
color6=[0.3010, 0.7450, 0.9330];
color7=[0.6350, 0.0780, 0.1840];

blue1=[0,0,128]/255;
blue2=[65,105,225]/255;
blue3=[135,206,250]/255;

red1=[178,34,34]/255;
red2=[220,20,60]/255;
red3=[250,128,114]/255;

L1=10000;
L2=20000;

% Rg 0.7 -----------------------------------------------------------------


N49071_Rg=importdata('N49_0.7_run1_Rg.txt');
N49072_Rg=importdata('N49_0.7_run2_Rg.txt');
N49073_Rg=importdata('N49_0.7_run3_Rg.txt');
N49074_Rg=importdata('N49_0.7_run4_Rg.txt');
N49075_Rg=importdata('N49_0.7_run5_Rg.txt');
N49076_Rg=importdata('N49_0.7_run6_Rg.txt');
N49077_Rg=importdata('N49_0.7_run7_Rg.txt');
N49078_Rg=importdata('N49_0.7_run8_Rg.txt');
N49079_Rg=importdata('N49_0.7_run9_Rg.txt');
N490710_Rg=importdata('N49_0.7_run10_Rg.txt');
N49071_Rg=N49071_Rg(L1:L2);
N49072_Rg=N49072_Rg(L1:L2);
N49073_Rg=N49073_Rg(L1:L2);
N49074_Rg=N49074_Rg(L1:L2);
N49075_Rg=N49075_Rg(L1:L2);
N49076_Rg=N49076_Rg(L1:L2);
N49077_Rg=N49077_Rg(L1:L2);
N49078_Rg=N49078_Rg(L1:L2);
N49079_Rg=N49079_Rg(L1:L2);
N490710_Rg=N490710_Rg(L1:L2);

Rg_07=mean([mean(N49071_Rg);mean(N49072_Rg);mean(N49073_Rg);mean(N49074_Rg);mean(N49075_Rg)]);
Rg_07_error=std([mean(N49071_Rg);mean(N49072_Rg);mean(N49073_Rg);mean(N49074_Rg);mean(N49075_Rg)]);

% Rg 0.75 -----------------------------------------------------------------

N490751_Rg=importdata('N49_0.75_run1_Rg.txt');
N490752_Rg=importdata('N49_0.75_run2_Rg.txt');
N490753_Rg=importdata('N49_0.75_run3_Rg.txt');
N490754_Rg=importdata('N49_0.75_run4_Rg.txt');
N490755_Rg=importdata('N49_0.75_run5_Rg.txt');
N490756_Rg=importdata('N49_0.75_run6_Rg.txt');
N490757_Rg=importdata('N49_0.75_run7_Rg.txt');
N490758_Rg=importdata('N49_0.75_run8_Rg.txt');
N490759_Rg=importdata('N49_0.75_run9_Rg.txt');
N4907510_Rg=importdata('N49_0.75_run10_Rg.txt');
N490751_Rg=N490751_Rg(L1:L2);
N490752_Rg=N490752_Rg(L1:L2);
N490753_Rg=N490753_Rg(L1:L2);
N490754_Rg=N490754_Rg(L1:L2);
N490755_Rg=N490755_Rg(L1:L2);
N490756_Rg=N490756_Rg(L1:L2);
N490757_Rg=N490757_Rg(L1:L2);
N490758_Rg=N490758_Rg(L1:L2);
N490759_Rg=N490759_Rg(L1:L2);
N4907510_Rg=N4907510_Rg(L1:L2);

Rg_075=mean([mean(N490751_Rg);mean(N490752_Rg);mean(N490753_Rg);mean(N490754_Rg);mean(N490755_Rg)]);
Rg_075_error=std([mean(N490751_Rg);mean(N490752_Rg);mean(N490753_Rg);mean(N490754_Rg);mean(N490755_Rg)]);

% Rg 0.8 -----------------------------------------------------------------

N49081_Rg=importdata('N49_0.8_run1_Rg.txt');
N49082_Rg=importdata('N49_0.8_run2_Rg.txt');
N49083_Rg=importdata('N49_0.8_run3_Rg.txt');
N49084_Rg=importdata('N49_0.8_run4_Rg.txt');
N49085_Rg=importdata('N49_0.8_run5_Rg.txt');
N49086_Rg=importdata('N49_0.8_run6_Rg.txt');
N49087_Rg=importdata('N49_0.8_run7_Rg.txt');
N49088_Rg=importdata('N49_0.8_run8_Rg.txt');
N49089_Rg=importdata('N49_0.8_run9_Rg.txt');
N490810_Rg=importdata('N49_0.8_run10_Rg.txt');
N49081_Rg=N49081_Rg(L1:L2);
N49082_Rg=N49082_Rg(L1:L2);
N49083_Rg=N49083_Rg(L1:L2);
N49084_Rg=N49084_Rg(L1:L2);
N49085_Rg=N49085_Rg(L1:L2);
N49086_Rg=N49086_Rg(L1:L2);
N49087_Rg=N49087_Rg(L1:L2);
N49088_Rg=N49088_Rg(L1:L2);
N49089_Rg=N49089_Rg(L1:L2);
N490810_Rg=N490810_Rg(L1:L2);

Rg_08=mean([mean(N49081_Rg);mean(N49082_Rg);mean(N49083_Rg);mean(N49084_Rg);mean(N49085_Rg)]);
Rg_08_error=std([mean(N49081_Rg);mean(N49082_Rg);mean(N49083_Rg);mean(N49084_Rg);mean(N49085_Rg)]);

% Rg 0.85 -----------------------------------------------------------------

N490851_Rg=importdata('N49_0.85_run1_Rg.txt');
N490852_Rg=importdata('N49_0.85_run2_Rg.txt');
N490853_Rg=importdata('N49_0.85_run3_Rg.txt');
N490854_Rg=importdata('N49_0.85_run4_Rg.txt');
N490855_Rg=importdata('N49_0.85_run5_Rg.txt');
N490856_Rg=importdata('N49_0.85_run6_Rg.txt');
N490857_Rg=importdata('N49_0.85_run7_Rg.txt');
N490858_Rg=importdata('N49_0.85_run8_Rg.txt');
N490859_Rg=importdata('N49_0.85_run9_Rg.txt');
N4908510_Rg=importdata('N49_0.85_run10_Rg.txt');
N490851_Rg=N490851_Rg(L1:L2);
N490852_Rg=N490852_Rg(L1:L2);
N490853_Rg=N490853_Rg(L1:L2);
N490854_Rg=N490854_Rg(L1:L2);
N490855_Rg=N490855_Rg(L1:L2);
N490856_Rg=N490856_Rg(L1:L2);
N490857_Rg=N490857_Rg(L1:L2);
N490858_Rg=N490858_Rg(L1:L2);
N490859_Rg=N490859_Rg(L1:L2);
N4908510_Rg=N4908510_Rg(L1:L2);

Rg_085=mean([mean(N490851_Rg);mean(N490852_Rg);mean(N490853_Rg);mean(N490854_Rg);mean(N490855_Rg)]);
Rg_085_error=std([mean(N490851_Rg);mean(N490852_Rg);mean(N490853_Rg);mean(N490854_Rg);mean(N490855_Rg)]);

% Rg 0.9 -----------------------------------------------------------------

N49091_Rg=importdata('N49_0.9_run1_Rg.txt');
N49092_Rg=importdata('N49_0.9_run2_Rg.txt');
N49093_Rg=importdata('N49_0.9_run3_Rg.txt');
N49094_Rg=importdata('N49_0.9_run4_Rg.txt');
N49095_Rg=importdata('N49_0.9_run5_Rg.txt');
N49096_Rg=importdata('N49_0.9_run6_Rg.txt');
N49097_Rg=importdata('N49_0.9_run7_Rg.txt');
N49098_Rg=importdata('N49_0.9_run8_Rg.txt');
N49099_Rg=importdata('N49_0.9_run9_Rg.txt');
N490910_Rg=importdata('N49_0.9_run10_Rg.txt');
N49091_Rg=N49091_Rg(L1:L2);
N49092_Rg=N49092_Rg(L1:L2);
N49093_Rg=N49093_Rg(L1:L2);
N49094_Rg=N49094_Rg(L1:L2);
N49095_Rg=N49095_Rg(L1:L2);
N49096_Rg=N49096_Rg(L1:L2);
N49097_Rg=N49097_Rg(L1:L2);
N49098_Rg=N49098_Rg(L1:L2);
N49099_Rg=N49099_Rg(L1:L2);
N490910_Rg=N490910_Rg(L1:L2);

Rg_09=mean([mean(N49091_Rg);mean(N49092_Rg);mean(N49093_Rg);mean(N49094_Rg);mean(N49095_Rg)]);
Rg_09_error=std([mean(N49091_Rg);mean(N49092_Rg);mean(N49093_Rg);mean(N49094_Rg);mean(N49095_Rg)]);

% Rg 0.95 -----------------------------------------------------------------

N490951_Rg=importdata('N49_0.95_run1_Rg.txt');
N490952_Rg=importdata('N49_0.95_run2_Rg.txt');
N490953_Rg=importdata('N49_0.95_run3_Rg.txt');
N490954_Rg=importdata('N49_0.95_run4_Rg.txt');
N490955_Rg=importdata('N49_0.95_run5_Rg.txt');
N490956_Rg=importdata('N49_0.95_run6_Rg.txt');
N490957_Rg=importdata('N49_0.95_run7_Rg.txt');
N490958_Rg=importdata('N49_0.95_run8_Rg.txt');
N490959_Rg=importdata('N49_0.95_run9_Rg.txt');
N4909510_Rg=importdata('N49_0.95_run10_Rg.txt');
N490951_Rg=N490951_Rg(L1:L2);
N490952_Rg=N490952_Rg(L1:L2);
N490953_Rg=N490953_Rg(L1:L2);
N490954_Rg=N490954_Rg(L1:L2);
N490955_Rg=N490955_Rg(L1:L2);
N490956_Rg=N490956_Rg(L1:L2);
N490957_Rg=N490957_Rg(L1:L2);
N490958_Rg=N490958_Rg(L1:L2);
N490959_Rg=N490959_Rg(L1:L2);
N4909510_Rg=N4909510_Rg(L1:L2);

Rg_095=mean([mean(N490951_Rg);mean(N490952_Rg);mean(N490953_Rg);mean(N490954_Rg);mean(N490955_Rg)]);
Rg_095_error=std([mean(N490951_Rg);mean(N490952_Rg);mean(N490953_Rg);mean(N490954_Rg);mean(N490955_Rg)]);

% Rg 1.0 -----------------------------------------------------------------

N49101_Rg=importdata('N49_1.0_run1_Rg.txt');
N49102_Rg=importdata('N49_1.0_run2_Rg.txt');
N49103_Rg=importdata('N49_1.0_run3_Rg.txt');
N49104_Rg=importdata('N49_1.0_run4_Rg.txt');
N49105_Rg=importdata('N49_1.0_run5_Rg.txt');
N49106_Rg=importdata('N49_1.0_run6_Rg.txt');
N49107_Rg=importdata('N49_1.0_run7_Rg.txt');
N49108_Rg=importdata('N49_1.0_run8_Rg.txt');
N49109_Rg=importdata('N49_1.0_run9_Rg.txt');
N491010_Rg=importdata('N49_1.0_run10_Rg.txt');
N49101_Rg=N49101_Rg(L1:L2);
N49102_Rg=N49102_Rg(L1:L2);
N49103_Rg=N49103_Rg(L1:L2);
N49104_Rg=N49104_Rg(L1:L2);
N49105_Rg=N49105_Rg(L1:L2);
N49106_Rg=N49106_Rg(L1:L2);
N49107_Rg=N49107_Rg(L1:L2);
N49108_Rg=N49108_Rg(L1:L2);
N49109_Rg=N49109_Rg(L1:L2);
N491010_Rg=N491010_Rg(L1:L2);

Rg_10=mean([mean(N49101_Rg);mean(N49102_Rg);mean(N49103_Rg);mean(N49104_Rg);mean(N49105_Rg)]);
Rg_10_error=std([mean(N49101_Rg);mean(N49102_Rg);mean(N49103_Rg);mean(N49104_Rg);mean(N49105_Rg)]);

% Rg 1.05 -----------------------------------------------------------------

N491051_Rg=importdata('N49_1.05_run1_Rg.txt');
N491052_Rg=importdata('N49_1.05_run2_Rg.txt');
N491053_Rg=importdata('N49_1.05_run3_Rg.txt');
N491054_Rg=importdata('N49_1.05_run4_Rg.txt');
N491055_Rg=importdata('N49_1.05_run5_Rg.txt');
N491056_Rg=importdata('N49_1.05_run6_Rg.txt');
N491057_Rg=importdata('N49_1.05_run7_Rg.txt');
N491058_Rg=importdata('N49_1.05_run8_Rg.txt');
N491059_Rg=importdata('N49_1.05_run9_Rg.txt');
N4910510_Rg=importdata('N49_1.05_run10_Rg.txt');
N491051_Rg=N491051_Rg(L1:L2);
N491052_Rg=N491052_Rg(L1:L2);
N491053_Rg=N491053_Rg(L1:L2);
N491054_Rg=N491054_Rg(L1:L2);
N491055_Rg=N491055_Rg(L1:L2);
N491056_Rg=N491056_Rg(L1:L2);
N491057_Rg=N491057_Rg(L1:L2);
N491058_Rg=N491058_Rg(L1:L2);
N491059_Rg=N491059_Rg(L1:L2);
N4910510_Rg=N4910510_Rg(L1:L2);

Rg_105=mean([mean(N491051_Rg);mean(N491052_Rg);mean(N491053_Rg);mean(N491054_Rg);mean(N491055_Rg)]);
Rg_105_error=std([mean(N491051_Rg);mean(N491052_Rg);mean(N491053_Rg);mean(N491054_Rg);mean(N491055_Rg)]);

% Rg 1.1 -----------------------------------------------------------------

N49111_Rg=importdata('N49_1.1_run1_Rg.txt');
N49112_Rg=importdata('N49_1.1_run2_Rg.txt');
N49113_Rg=importdata('N49_1.1_run3_Rg.txt');
N49114_Rg=importdata('N49_1.1_run4_Rg.txt');
N49115_Rg=importdata('N49_1.1_run5_Rg.txt');
N49116_Rg=importdata('N49_1.1_run6_Rg.txt');
N49117_Rg=importdata('N49_1.1_run7_Rg.txt');
N49118_Rg=importdata('N49_1.1_run8_Rg.txt');
N49119_Rg=importdata('N49_1.1_run9_Rg.txt');
N491110_Rg=importdata('N49_1.1_run10_Rg.txt');
N49111_Rg=N49111_Rg(L1:L2);
N49112_Rg=N49112_Rg(L1:L2);
N49113_Rg=N49113_Rg(L1:L2);
N49114_Rg=N49114_Rg(L1:L2);
N49115_Rg=N49115_Rg(L1:L2);
N49116_Rg=N49116_Rg(L1:L2);
N49117_Rg=N49117_Rg(L1:L2);
N49118_Rg=N49118_Rg(L1:L2);
N49119_Rg=N49119_Rg(L1:L2);
N491110_Rg=N491110_Rg(L1:L2);

Rg_11=mean([mean(N49111_Rg);mean(N49112_Rg);mean(N49113_Rg);mean(N49114_Rg);mean(N49115_Rg)]);
Rg_11_error=std([mean(N49111_Rg);mean(N49112_Rg);mean(N49113_Rg);mean(N49114_Rg);mean(N49115_Rg)]);

% Plot data ----------------------------------------------------------------------------------------------------------------------

x=[1.0;20];
y=[1.0;20];

% calculate the shift as a funciton of eps for both Rg and FRET
shift_Rg=zeros(9,2);
shift_Rg(:,1)=[0.7:0.05:1.1];
shift_Rg(1,2)=Rg_07-Rg_105;
shift_Rg(2,2)=Rg_075-Rg_105;
shift_Rg(3,2)=Rg_08-Rg_105;
shift_Rg(4,2)=Rg_085-Rg_105;
shift_Rg(5,2)=Rg_09-Rg_105;
shift_Rg(6,2)=Rg_095-Rg_105;
shift_Rg(7,2)=Rg_10-Rg_105;
shift_Rg(8,2)=Rg_105-Rg_105;
shift_Rg(9,2)=Rg_11-Rg_105;

Rg_N49_shift=[0.1078;0.1078];
Rg_N49_shift2=[0.1593;0.1593];

x=[0.7;1.1];

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(x,Rg_N49_shift,'--','Linewidth',3,'Color',color2,'MarkerSize',12,'MarkerFaceColor',color2);
plot(x,Rg_N49_shift2,'--','Linewidth',3,'Color',color5,'MarkerSize',12,'MarkerFaceColor',color5);
plot(shift_Rg(:,1),shift_Rg(:,2),'o-','Linewidth',3,'Color',color1,'MarkerSize',12,'MarkerFaceColor',color1);
legend({' 1.5%',' 5.0%'},'location','southwest');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([0.7 1.1 -0.05 0.2]);
legend boxoff;
box on;
print(fig,'Rg_shift.png','-dpng');

N49_Rg_error=abs( shift_Rg(:,2)-Rg_N49_shift(1) );
N49_Rg_error2=abs( shift_Rg(:,2)-Rg_N49_shift2(1) );


