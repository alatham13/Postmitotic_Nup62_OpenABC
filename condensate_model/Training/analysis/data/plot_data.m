
clear all;
close all;

color1=[0, 0.4470, 0.7410];
color2=[0.8500, 0.3250, 0.0980];
color3=[0.9290, 0.6940, 0.1250];
color4=[0.4940, 0.1840, 0.5560];
color5=[0.4660, 0.6740, 0.1880];
color6=[0.3010, 0.7450, 0.9330];
color7=[0.6350, 0.0780, 0.1840];

% Experimental Rg from the unlabeled protein
% Order: N49, N98, NSP, NUL, NUS
exp_Rg=[1.59;2.86;4.1;3.0;2.49];
exp_Rg_error=[0.13;0.13;0.3;0.3;0.13];

% Experimental FRET from the labeled protein
% Order: N49, N98, NSP, NUL, NUS
exp_FRET=[3.6;5.3;6.9;6.6;6.2];
exp_FRET_error=zeros(5,1)+0.3;
N=5;
L1=10000;
L2=20000;

% Rg 0.9 -----------------------------------------------------------------

Rg_09=zeros(5,1);
Rg_09_error=zeros(5,1);

N49091_Rg=importdata('N49_0.9_run1_Rg.txt');
N49092_Rg=importdata('N49_0.9_run2_Rg.txt');
N49093_Rg=importdata('N49_0.9_run3_Rg.txt');
N49094_Rg=importdata('N49_0.9_run4_Rg.txt');
N49095_Rg=importdata('N49_0.9_run5_Rg.txt');
N49091_Rg=N49091_Rg(L1:L2);
N49092_Rg=N49092_Rg(L1:L2);
N49093_Rg=N49093_Rg(L1:L2);
N49094_Rg=N49094_Rg(L1:L2);
N49095_Rg=N49095_Rg(L1:L2);

Rg_09(1)=mean([mean(N49091_Rg);mean(N49092_Rg);mean(N49093_Rg);mean(N49094_Rg);mean(N49095_Rg)]);
Rg_09_error(1)=std([mean(N49091_Rg);mean(N49092_Rg);mean(N49093_Rg);mean(N49094_Rg);mean(N49095_Rg)]);


N98091_Rg=importdata('N98_0.9_run1_Rg.txt');
N98092_Rg=importdata('N98_0.9_run2_Rg.txt');
N98093_Rg=importdata('N98_0.9_run3_Rg.txt');
N98094_Rg=importdata('N98_0.9_run4_Rg.txt');
N98095_Rg=importdata('N98_0.9_run5_Rg.txt');
N98091_Rg=N98091_Rg(L1:L2);
N98092_Rg=N98092_Rg(L1:L2);
N98093_Rg=N98093_Rg(L1:L2);
N98094_Rg=N98094_Rg(L1:L2);
N98095_Rg=N98095_Rg(L1:L2);

Rg_09(2)=mean([mean(N98091_Rg);mean(N98092_Rg);mean(N98093_Rg);mean(N98094_Rg);mean(N98095_Rg)]);
Rg_09_error(2)=std([mean(N98091_Rg);mean(N98092_Rg);mean(N98093_Rg);mean(N98094_Rg);mean(N98095_Rg)]);

NSP091_Rg=importdata('NSP_0.9_run1_Rg.txt');
NSP092_Rg=importdata('NSP_0.9_run2_Rg.txt');
NSP093_Rg=importdata('NSP_0.9_run3_Rg.txt');
NSP094_Rg=importdata('NSP_0.9_run4_Rg.txt');
NSP095_Rg=importdata('NSP_0.9_run5_Rg.txt');
NSP091_Rg=NSP091_Rg(L1:L2);
NSP092_Rg=NSP092_Rg(L1:L2);
NSP093_Rg=NSP093_Rg(L1:L2);
NSP094_Rg=NSP094_Rg(L1:L2);
NSP095_Rg=NSP095_Rg(L1:L2);

Rg_09(3)=mean([mean(NSP091_Rg);mean(NSP092_Rg);mean(NSP093_Rg);mean(NSP094_Rg);mean(NSP095_Rg)]);
Rg_09_error(3)=std([mean(NSP091_Rg);mean(NSP092_Rg);mean(NSP093_Rg);mean(NSP094_Rg);mean(NSP095_Rg)]);

NUL091_Rg=importdata('NUL_0.9_run1_Rg.txt');
NUL092_Rg=importdata('NUL_0.9_run2_Rg.txt');
NUL093_Rg=importdata('NUL_0.9_run3_Rg.txt');
NUL094_Rg=importdata('NUL_0.9_run4_Rg.txt');
NUL095_Rg=importdata('NUL_0.9_run5_Rg.txt');
NUL091_Rg=NUL091_Rg(L1:L2);
NUL092_Rg=NUL092_Rg(L1:L2);
NUL093_Rg=NUL093_Rg(L1:L2);
NUL094_Rg=NUL094_Rg(L1:L2);
NUL095_Rg=NUL095_Rg(L1:L2);

Rg_09(4)=mean([mean(NUL091_Rg);mean(NUL092_Rg);mean(NUL093_Rg);mean(NUL094_Rg);mean(NUL095_Rg)]);
Rg_09_error(4)=std([mean(NUL091_Rg);mean(NUL092_Rg);mean(NUL093_Rg);mean(NUL094_Rg);mean(NUL095_Rg)]);

NUS091_Rg=importdata('NUS_0.9_run1_Rg.txt');
NUS092_Rg=importdata('NUS_0.9_run2_Rg.txt');
NUS093_Rg=importdata('NUS_0.9_run3_Rg.txt');
NUS094_Rg=importdata('NUS_0.9_run4_Rg.txt');
NUS095_Rg=importdata('NUS_0.9_run5_Rg.txt');
NUS091_Rg=NUS091_Rg(L1:L2);
NUS092_Rg=NUS092_Rg(L1:L2);
NUS093_Rg=NUS093_Rg(L1:L2);
NUS094_Rg=NUS094_Rg(L1:L2);
NUS095_Rg=NUS095_Rg(L1:L2);

Rg_09(5)=mean([mean(NUS091_Rg);mean(NUS092_Rg);mean(NUS093_Rg);mean(NUS094_Rg);mean(NUS095_Rg)]);
Rg_09_error(5)=std([mean(NUS091_Rg);mean(NUS092_Rg);mean(NUS093_Rg);mean(NUS094_Rg);mean(NUS095_Rg)]);


% FRET 0.9 -----------------------------------------------------------------

FRET_09=zeros(5,1);
FRET_09_error=zeros(5,1);

N49091_FRET=importdata('N49_0.9_run1_FRET.txt');
N49092_FRET=importdata('N49_0.9_run2_FRET.txt');
N49093_FRET=importdata('N49_0.9_run3_FRET.txt');
N49094_FRET=importdata('N49_0.9_run4_FRET.txt');
N49095_FRET=importdata('N49_0.9_run5_FRET.txt');
N49091_FRET=N49091_FRET(L1:L2);
N49092_FRET=N49092_FRET(L1:L2);
N49093_FRET=N49093_FRET(L1:L2);
N49094_FRET=N49094_FRET(L1:L2);
N49095_FRET=N49095_FRET(L1:L2);

FRET_09(1)=mean([mean(N49091_FRET);mean(N49092_FRET);mean(N49093_FRET);mean(N49094_FRET);mean(N49095_FRET)]);
FRET_09_error(1)=std([mean(N49091_FRET);mean(N49092_FRET);mean(N49093_FRET);mean(N49094_FRET);mean(N49095_FRET)]);


N98091_FRET=importdata('N98_0.9_run1_FRET.txt');
N98092_FRET=importdata('N98_0.9_run2_FRET.txt');
N98093_FRET=importdata('N98_0.9_run3_FRET.txt');
N98094_FRET=importdata('N98_0.9_run4_FRET.txt');
N98095_FRET=importdata('N98_0.9_run5_FRET.txt');
N98091_FRET=N98091_FRET(L1:L2);
N98092_FRET=N98092_FRET(L1:L2);
N98093_FRET=N98093_FRET(L1:L2);
N98094_FRET=N98094_FRET(L1:L2);
N98095_FRET=N98095_FRET(L1:L2);

FRET_09(2)=mean([mean(N98091_FRET);mean(N98092_FRET);mean(N98093_FRET);mean(N98094_FRET);mean(N98095_FRET)]);
FRET_09_error(2)=std([mean(N98091_FRET);mean(N98092_FRET);mean(N98093_FRET);mean(N98094_FRET);mean(N98095_FRET)]);

NSP091_FRET=importdata('NSP_0.9_run1_FRET.txt');
NSP092_FRET=importdata('NSP_0.9_run2_FRET.txt');
NSP093_FRET=importdata('NSP_0.9_run3_FRET.txt');
NSP094_FRET=importdata('NSP_0.9_run4_FRET.txt');
NSP095_FRET=importdata('NSP_0.9_run5_FRET.txt');
NSP091_FRET=NSP091_FRET(L1:L2);
NSP092_FRET=NSP092_FRET(L1:L2);
NSP093_FRET=NSP093_FRET(L1:L2);
NSP094_FRET=NSP094_FRET(L1:L2);
NSP095_FRET=NSP095_FRET(L1:L2);

FRET_09(3)=mean([mean(NSP091_FRET);mean(NSP092_FRET);mean(NSP093_FRET);mean(NSP094_FRET);mean(NSP095_FRET)]);
FRET_09_error(3)=std([mean(NSP091_FRET);mean(NSP092_FRET);mean(NSP093_FRET);mean(NSP094_FRET);mean(NSP095_FRET)]);

NUL091_FRET=importdata('NUL_0.9_run1_FRET.txt');
NUL092_FRET=importdata('NUL_0.9_run2_FRET.txt');
NUL093_FRET=importdata('NUL_0.9_run3_FRET.txt');
NUL094_FRET=importdata('NUL_0.9_run4_FRET.txt');
NUL095_FRET=importdata('NUL_0.9_run5_FRET.txt');
NUL091_FRET=NUL091_FRET(L1:L2);
NUL092_FRET=NUL092_FRET(L1:L2);
NUL093_FRET=NUL093_FRET(L1:L2);
NUL094_FRET=NUL094_FRET(L1:L2);
NUL095_FRET=NUL095_FRET(L1:L2);

FRET_09(4)=mean([mean(NUL091_FRET);mean(NUL092_FRET);mean(NUL093_FRET);mean(NUL094_FRET);mean(NUL095_FRET)]);
FRET_09_error(4)=std([mean(NUL091_FRET);mean(NUL092_FRET);mean(NUL093_FRET);mean(NUL094_FRET);mean(NUL095_FRET)]);

NUS091_FRET=importdata('NUS_0.9_run1_FRET.txt');
NUS092_FRET=importdata('NUS_0.9_run2_FRET.txt');
NUS093_FRET=importdata('NUS_0.9_run3_FRET.txt');
NUS094_FRET=importdata('NUS_0.9_run4_FRET.txt');
NUS095_FRET=importdata('NUS_0.9_run5_FRET.txt');
NUS091_FRET=NUS091_FRET(L1:L2);
NUS092_FRET=NUS092_FRET(L1:L2);
NUS093_FRET=NUS093_FRET(L1:L2);
NUS094_FRET=NUS094_FRET(L1:L2);
NUS095_FRET=NUS095_FRET(L1:L2);

FRET_09(5)=mean([mean(NUS091_FRET);mean(NUS092_FRET);mean(NUS093_FRET);mean(NUS094_FRET);mean(NUS095_FRET)]);
FRET_09_error(5)=std([mean(NUS091_FRET);mean(NUS092_FRET);mean(NUS093_FRET);mean(NUS094_FRET);mean(NUS095_FRET)]);

% Rg 0.95 -----------------------------------------------------------------

Rg_095=zeros(5,1);
Rg_095_error=zeros(5,1);

N490951_Rg=importdata('N49_0.95_run1_Rg.txt');
N490952_Rg=importdata('N49_0.95_run2_Rg.txt');
N490953_Rg=importdata('N49_0.95_run3_Rg.txt');
N490954_Rg=importdata('N49_0.95_run4_Rg.txt');
N490955_Rg=importdata('N49_0.95_run5_Rg.txt');
N490951_Rg=N490951_Rg(L1:L2);
N490952_Rg=N490952_Rg(L1:L2);
N490953_Rg=N490953_Rg(L1:L2);
N490954_Rg=N490954_Rg(L1:L2);
N490955_Rg=N490955_Rg(L1:L2);

Rg_095(1)=mean([mean(N490951_Rg);mean(N490952_Rg);mean(N490953_Rg);mean(N490954_Rg);mean(N490955_Rg)]);
Rg_095_error(1)=std([mean(N490951_Rg);mean(N490952_Rg);mean(N490953_Rg);mean(N490954_Rg);mean(N490955_Rg)]);


N980951_Rg=importdata('N98_0.95_run1_Rg.txt');
N980952_Rg=importdata('N98_0.95_run2_Rg.txt');
N980953_Rg=importdata('N98_0.95_run3_Rg.txt');
N980954_Rg=importdata('N98_0.95_run4_Rg.txt');
N980955_Rg=importdata('N98_0.95_run5_Rg.txt');
N980951_Rg=N980951_Rg(L1:L2);
N980952_Rg=N980952_Rg(L1:L2);
N980953_Rg=N980953_Rg(L1:L2);
N980954_Rg=N980954_Rg(L1:L2);
N980955_Rg=N980955_Rg(L1:L2);

Rg_095(2)=mean([mean(N980951_Rg);mean(N980952_Rg);mean(N980953_Rg);mean(N980954_Rg);mean(N980955_Rg)]);
Rg_095_error(2)=std([mean(N980951_Rg);mean(N980952_Rg);mean(N980953_Rg);mean(N980954_Rg);mean(N980955_Rg)]);

NSP0951_Rg=importdata('NSP_0.95_run1_Rg.txt');
NSP0952_Rg=importdata('NSP_0.95_run2_Rg.txt');
NSP0953_Rg=importdata('NSP_0.95_run3_Rg.txt');
NSP0954_Rg=importdata('NSP_0.95_run4_Rg.txt');
NSP0955_Rg=importdata('NSP_0.95_run5_Rg.txt');
NSP0951_Rg=NSP0951_Rg(L1:L2);
NSP0952_Rg=NSP0952_Rg(L1:L2);
NSP0953_Rg=NSP0953_Rg(L1:L2);
NSP0954_Rg=NSP0954_Rg(L1:L2);
NSP0955_Rg=NSP0955_Rg(L1:L2);

Rg_095(3)=mean([mean(NSP0951_Rg);mean(NSP0952_Rg);mean(NSP0953_Rg);mean(NSP0954_Rg);mean(NSP0955_Rg)]);
Rg_095_error(3)=std([mean(NSP0951_Rg);mean(NSP0952_Rg);mean(NSP0953_Rg);mean(NSP0954_Rg);mean(NSP0955_Rg)]);

NUL0951_Rg=importdata('NUL_0.95_run1_Rg.txt');
NUL0952_Rg=importdata('NUL_0.95_run2_Rg.txt');
NUL0953_Rg=importdata('NUL_0.95_run3_Rg.txt');
NUL0954_Rg=importdata('NUL_0.95_run4_Rg.txt');
NUL0955_Rg=importdata('NUL_0.95_run5_Rg.txt');
NUL0951_Rg=NUL0951_Rg(L1:L2);
NUL0952_Rg=NUL0952_Rg(L1:L2);
NUL0953_Rg=NUL0953_Rg(L1:L2);
NUL0954_Rg=NUL0954_Rg(L1:L2);
NUL0955_Rg=NUL0955_Rg(L1:L2);

Rg_095(4)=mean([mean(NUL0951_Rg);mean(NUL0952_Rg);mean(NUL0953_Rg);mean(NUL0954_Rg);mean(NUL0955_Rg)]);
Rg_095_error(4)=std([mean(NUL0951_Rg);mean(NUL0952_Rg);mean(NUL0953_Rg);mean(NUL0954_Rg);mean(NUL0955_Rg)]);

NUS0951_Rg=importdata('NUS_0.95_run1_Rg.txt');
NUS0952_Rg=importdata('NUS_0.95_run2_Rg.txt');
NUS0953_Rg=importdata('NUS_0.95_run3_Rg.txt');
NUS0954_Rg=importdata('NUS_0.95_run4_Rg.txt');
NUS0955_Rg=importdata('NUS_0.95_run5_Rg.txt');
NUS0951_Rg=NUS0951_Rg(L1:L2);
NUS0952_Rg=NUS0952_Rg(L1:L2);
NUS0953_Rg=NUS0953_Rg(L1:L2);
NUS0954_Rg=NUS0954_Rg(L1:L2);
NUS0955_Rg=NUS0955_Rg(L1:L2);

Rg_095(5)=mean([mean(NUS0951_Rg);mean(NUS0952_Rg);mean(NUS0953_Rg);mean(NUS0954_Rg);mean(NUS0955_Rg)]);
Rg_095_error(5)=std([mean(NUS0951_Rg);mean(NUS0952_Rg);mean(NUS0953_Rg);mean(NUS0954_Rg);mean(NUS0955_Rg)]);


% FRET 0.95 -----------------------------------------------------------------

FRET_095=zeros(5,1);
FRET_095_error=zeros(5,1);

N490951_FRET=importdata('N49_0.95_run1_FRET.txt');
N490952_FRET=importdata('N49_0.95_run2_FRET.txt');
N490953_FRET=importdata('N49_0.95_run3_FRET.txt');
N490954_FRET=importdata('N49_0.95_run4_FRET.txt');
N490955_FRET=importdata('N49_0.95_run5_FRET.txt');
N490951_FRET=N490951_FRET(L1:L2);
N490952_FRET=N490952_FRET(L1:L2);
N490953_FRET=N490953_FRET(L1:L2);
N490954_FRET=N490954_FRET(L1:L2);
N490955_FRET=N490955_FRET(L1:L2);

FRET_095(1)=mean([mean(N490951_FRET);mean(N490952_FRET);mean(N490953_FRET);mean(N490954_FRET);mean(N490955_FRET)]);
FRET_095_error(1)=std([mean(N490951_FRET);mean(N490952_FRET);mean(N490953_FRET);mean(N490954_FRET);mean(N490955_FRET)]);


N980951_FRET=importdata('N98_0.95_run1_FRET.txt');
N980952_FRET=importdata('N98_0.95_run2_FRET.txt');
N980953_FRET=importdata('N98_0.95_run3_FRET.txt');
N980954_FRET=importdata('N98_0.95_run4_FRET.txt');
N980955_FRET=importdata('N98_0.95_run5_FRET.txt');
N980951_FRET=N980951_FRET(L1:L2);
N980952_FRET=N980952_FRET(L1:L2);
N980953_FRET=N980953_FRET(L1:L2);
N980954_FRET=N980954_FRET(L1:L2);
N980955_FRET=N980955_FRET(L1:L2);

FRET_095(2)=mean([mean(N980951_FRET);mean(N980952_FRET);mean(N980953_FRET);mean(N980954_FRET);mean(N980955_FRET)]);
FRET_095_error(2)=std([mean(N980951_FRET);mean(N980952_FRET);mean(N980953_FRET);mean(N980954_FRET);mean(N980955_FRET)]);

NSP0951_FRET=importdata('NSP_0.95_run1_FRET.txt');
NSP0952_FRET=importdata('NSP_0.95_run2_FRET.txt');
NSP0953_FRET=importdata('NSP_0.95_run3_FRET.txt');
NSP0954_FRET=importdata('NSP_0.95_run4_FRET.txt');
NSP0955_FRET=importdata('NSP_0.95_run5_FRET.txt');
NSP0951_FRET=NSP0951_FRET(L1:L2);
NSP0952_FRET=NSP0952_FRET(L1:L2);
NSP0953_FRET=NSP0953_FRET(L1:L2);
NSP0954_FRET=NSP0954_FRET(L1:L2);
NSP0955_FRET=NSP0955_FRET(L1:L2);

FRET_095(3)=mean([mean(NSP0951_FRET);mean(NSP0952_FRET);mean(NSP0953_FRET);mean(NSP0954_FRET);mean(NSP0955_FRET)]);
FRET_095_error(3)=std([mean(NSP0951_FRET);mean(NSP0952_FRET);mean(NSP0953_FRET);mean(NSP0954_FRET);mean(NSP0955_FRET)]);

NUL0951_FRET=importdata('NUL_0.95_run1_FRET.txt');
NUL0952_FRET=importdata('NUL_0.95_run2_FRET.txt');
NUL0953_FRET=importdata('NUL_0.95_run3_FRET.txt');
NUL0954_FRET=importdata('NUL_0.95_run4_FRET.txt');
NUL0955_FRET=importdata('NUL_0.95_run5_FRET.txt');
NUL0951_FRET=NUL0951_FRET(L1:L2);
NUL0952_FRET=NUL0952_FRET(L1:L2);
NUL0953_FRET=NUL0953_FRET(L1:L2);
NUL0954_FRET=NUL0954_FRET(L1:L2);
NUL0955_FRET=NUL0955_FRET(L1:L2);

FRET_095(4)=mean([mean(NUL0951_FRET);mean(NUL0952_FRET);mean(NUL0953_FRET);mean(NUL0954_FRET);mean(NUL0955_FRET)]);
FRET_095_error(4)=std([mean(NUL0951_FRET);mean(NUL0952_FRET);mean(NUL0953_FRET);mean(NUL0954_FRET);mean(NUL0955_FRET)]);

NUS0951_FRET=importdata('NUS_0.95_run1_FRET.txt');
NUS0952_FRET=importdata('NUS_0.95_run2_FRET.txt');
NUS0953_FRET=importdata('NUS_0.95_run3_FRET.txt');
NUS0954_FRET=importdata('NUS_0.95_run4_FRET.txt');
NUS0955_FRET=importdata('NUS_0.95_run5_FRET.txt');
NUS0951_FRET=NUS0951_FRET(L1:L2);
NUS0952_FRET=NUS0952_FRET(L1:L2);
NUS0953_FRET=NUS0953_FRET(L1:L2);
NUS0954_FRET=NUS0954_FRET(L1:L2);
NUS0955_FRET=NUS0955_FRET(L1:L2);

FRET_095(5)=mean([mean(NUS0951_FRET);mean(NUS0952_FRET);mean(NUS0953_FRET);mean(NUS0954_FRET);mean(NUS0955_FRET)]);
FRET_095_error(5)=std([mean(NUS0951_FRET);mean(NUS0952_FRET);mean(NUS0953_FRET);mean(NUS0954_FRET);mean(NUS0955_FRET)]);

% Rg 1.0 -----------------------------------------------------------------

Rg_10=zeros(5,1);
Rg_10_error=zeros(5,1);

N49101_Rg=importdata('N49_1.0_run1_Rg.txt');
N49102_Rg=importdata('N49_1.0_run2_Rg.txt');
N49103_Rg=importdata('N49_1.0_run3_Rg.txt');
N49104_Rg=importdata('N49_1.0_run4_Rg.txt');
N49105_Rg=importdata('N49_1.0_run5_Rg.txt');
N49101_Rg=N49101_Rg(L1:L2);
N49102_Rg=N49102_Rg(L1:L2);
N49103_Rg=N49103_Rg(L1:L2);
N49104_Rg=N49104_Rg(L1:L2);
N49105_Rg=N49105_Rg(L1:L2);

Rg_10(1)=mean([mean(N49101_Rg);mean(N49102_Rg);mean(N49103_Rg);mean(N49104_Rg);mean(N49105_Rg)]);
Rg_10_error(1)=std([mean(N49101_Rg);mean(N49102_Rg);mean(N49103_Rg);mean(N49104_Rg);mean(N49105_Rg)]);


N98101_Rg=importdata('N98_1.0_run1_Rg.txt');
N98102_Rg=importdata('N98_1.0_run2_Rg.txt');
N98103_Rg=importdata('N98_1.0_run3_Rg.txt');
N98104_Rg=importdata('N98_1.0_run4_Rg.txt');
N98105_Rg=importdata('N98_1.0_run5_Rg.txt');
N98101_Rg=N98101_Rg(L1:L2);
N98102_Rg=N98102_Rg(L1:L2);
N98103_Rg=N98103_Rg(L1:L2);
N98104_Rg=N98104_Rg(L1:L2);
N98105_Rg=N98105_Rg(L1:L2);

Rg_10(2)=mean([mean(N98101_Rg);mean(N98102_Rg);mean(N98103_Rg);mean(N98104_Rg);mean(N98105_Rg)]);
Rg_10_error(2)=std([mean(N98101_Rg);mean(N98102_Rg);mean(N98103_Rg);mean(N98104_Rg);mean(N98105_Rg)]);

NSP101_Rg=importdata('NSP_1.0_run1_Rg.txt');
NSP102_Rg=importdata('NSP_1.0_run2_Rg.txt');
NSP103_Rg=importdata('NSP_1.0_run3_Rg.txt');
NSP104_Rg=importdata('NSP_1.0_run4_Rg.txt');
NSP105_Rg=importdata('NSP_1.0_run5_Rg.txt');
NSP101_Rg=NSP101_Rg(L1:L2);
NSP102_Rg=NSP102_Rg(L1:L2);
NSP103_Rg=NSP103_Rg(L1:L2);
NSP104_Rg=NSP104_Rg(L1:L2);
NSP105_Rg=NSP105_Rg(L1:L2);

Rg_10(3)=mean([mean(NSP101_Rg);mean(NSP102_Rg);mean(NSP103_Rg);mean(NSP104_Rg);mean(NSP105_Rg)]);
Rg_10_error(3)=std([mean(NSP101_Rg);mean(NSP102_Rg);mean(NSP103_Rg);mean(NSP104_Rg);mean(NSP105_Rg)]);

NUL101_Rg=importdata('NUL_1.0_run1_Rg.txt');
NUL102_Rg=importdata('NUL_1.0_run2_Rg.txt');
NUL103_Rg=importdata('NUL_1.0_run3_Rg.txt');
NUL104_Rg=importdata('NUL_1.0_run4_Rg.txt');
NUL105_Rg=importdata('NUL_1.0_run5_Rg.txt');
NUL101_Rg=NUL101_Rg(L1:L2);
NUL102_Rg=NUL102_Rg(L1:L2);
NUL103_Rg=NUL103_Rg(L1:L2);
NUL104_Rg=NUL104_Rg(L1:L2);
NUL105_Rg=NUL105_Rg(L1:L2);

Rg_10(4)=mean([mean(NUL101_Rg);mean(NUL102_Rg);mean(NUL103_Rg);mean(NUL104_Rg);mean(NUL105_Rg)]);
Rg_10_error(4)=std([mean(NUL101_Rg);mean(NUL102_Rg);mean(NUL103_Rg);mean(NUL104_Rg);mean(NUL105_Rg)]);

NUS101_Rg=importdata('NUS_1.0_run1_Rg.txt');
NUS102_Rg=importdata('NUS_1.0_run2_Rg.txt');
NUS103_Rg=importdata('NUS_1.0_run3_Rg.txt');
NUS104_Rg=importdata('NUS_1.0_run4_Rg.txt');
NUS105_Rg=importdata('NUS_1.0_run5_Rg.txt');
NUS101_Rg=NUS101_Rg(L1:L2);
NUS102_Rg=NUS102_Rg(L1:L2);
NUS103_Rg=NUS103_Rg(L1:L2);
NUS104_Rg=NUS104_Rg(L1:L2);
NUS105_Rg=NUS105_Rg(L1:L2);

Rg_10(5)=mean([mean(NUS101_Rg);mean(NUS102_Rg);mean(NUS103_Rg);mean(NUS104_Rg);mean(NUS105_Rg)]);
Rg_10_error(5)=std([mean(NUS101_Rg);mean(NUS102_Rg);mean(NUS103_Rg);mean(NUS104_Rg);mean(NUS105_Rg)]);


% FRET 1.0 -----------------------------------------------------------------

FRET_10=zeros(5,1);
FRET_10_error=zeros(5,1);

N49101_FRET=importdata('N49_1.0_run1_FRET.txt');
N49102_FRET=importdata('N49_1.0_run2_FRET.txt');
N49103_FRET=importdata('N49_1.0_run3_FRET.txt');
N49104_FRET=importdata('N49_1.0_run4_FRET.txt');
N49105_FRET=importdata('N49_1.0_run5_FRET.txt');
N49101_FRET=N49101_FRET(L1:L2);
N49102_FRET=N49102_FRET(L1:L2);
N49103_FRET=N49103_FRET(L1:L2);
N49104_FRET=N49104_FRET(L1:L2);
N49105_FRET=N49105_FRET(L1:L2);

FRET_10(1)=mean([mean(N49101_FRET);mean(N49102_FRET);mean(N49103_FRET);mean(N49104_FRET);mean(N49105_FRET)]);
FRET_10_error(1)=std([mean(N49101_FRET);mean(N49102_FRET);mean(N49103_FRET);mean(N49104_FRET);mean(N49105_FRET)]);


N98101_FRET=importdata('N98_1.0_run1_FRET.txt');
N98102_FRET=importdata('N98_1.0_run2_FRET.txt');
N98103_FRET=importdata('N98_1.0_run3_FRET.txt');
N98104_FRET=importdata('N98_1.0_run4_FRET.txt');
N98105_FRET=importdata('N98_1.0_run5_FRET.txt');
N98101_FRET=N98101_FRET(L1:L2);
N98102_FRET=N98102_FRET(L1:L2);
N98103_FRET=N98103_FRET(L1:L2);
N98104_FRET=N98104_FRET(L1:L2);
N98105_FRET=N98105_FRET(L1:L2);

FRET_10(2)=mean([mean(N98101_FRET);mean(N98102_FRET);mean(N98103_FRET);mean(N98104_FRET);mean(N98105_FRET)]);
FRET_10_error(2)=std([mean(N98101_FRET);mean(N98102_FRET);mean(N98103_FRET);mean(N98104_FRET);mean(N98105_FRET)]);

NSP101_FRET=importdata('NSP_1.0_run1_FRET.txt');
NSP102_FRET=importdata('NSP_1.0_run2_FRET.txt');
NSP103_FRET=importdata('NSP_1.0_run3_FRET.txt');
NSP104_FRET=importdata('NSP_1.0_run4_FRET.txt');
NSP105_FRET=importdata('NSP_1.0_run5_FRET.txt');
NSP101_FRET=NSP101_FRET(L1:L2);
NSP102_FRET=NSP102_FRET(L1:L2);
NSP103_FRET=NSP103_FRET(L1:L2);
NSP104_FRET=NSP104_FRET(L1:L2);
NSP105_FRET=NSP105_FRET(L1:L2);

FRET_10(3)=mean([mean(NSP101_FRET);mean(NSP102_FRET);mean(NSP103_FRET);mean(NSP104_FRET);mean(NSP105_FRET)]);
FRET_10_error(3)=std([mean(NSP101_FRET);mean(NSP102_FRET);mean(NSP103_FRET);mean(NSP104_FRET);mean(NSP105_FRET)]);

NUL101_FRET=importdata('NUL_1.0_run1_FRET.txt');
NUL102_FRET=importdata('NUL_1.0_run2_FRET.txt');
NUL103_FRET=importdata('NUL_1.0_run3_FRET.txt');
NUL104_FRET=importdata('NUL_1.0_run4_FRET.txt');
NUL105_FRET=importdata('NUL_1.0_run5_FRET.txt');
NUL101_FRET=NUL101_FRET(L1:L2);
NUL102_FRET=NUL102_FRET(L1:L2);
NUL103_FRET=NUL103_FRET(L1:L2);
NUL104_FRET=NUL104_FRET(L1:L2);
NUL105_FRET=NUL105_FRET(L1:L2);

FRET_10(4)=mean([mean(NUL101_FRET);mean(NUL102_FRET);mean(NUL103_FRET);mean(NUL104_FRET);mean(NUL105_FRET)]);
FRET_10_error(4)=std([mean(NUL101_FRET);mean(NUL102_FRET);mean(NUL103_FRET);mean(NUL104_FRET);mean(NUL105_FRET)]);

NUS101_FRET=importdata('NUS_1.0_run1_FRET.txt');
NUS102_FRET=importdata('NUS_1.0_run2_FRET.txt');
NUS103_FRET=importdata('NUS_1.0_run3_FRET.txt');
NUS104_FRET=importdata('NUS_1.0_run4_FRET.txt');
NUS105_FRET=importdata('NUS_1.0_run5_FRET.txt');
NUS101_FRET=NUS101_FRET(L1:L2);
NUS102_FRET=NUS102_FRET(L1:L2);
NUS103_FRET=NUS103_FRET(L1:L2);
NUS104_FRET=NUS104_FRET(L1:L2);
NUS105_FRET=NUS105_FRET(L1:L2);

FRET_10(5)=mean([mean(NUS101_FRET);mean(NUS102_FRET);mean(NUS103_FRET);mean(NUS104_FRET);mean(NUS105_FRET)]);
FRET_10_error(5)=std([mean(NUS101_FRET);mean(NUS102_FRET);mean(NUS103_FRET);mean(NUS104_FRET);mean(NUS105_FRET)]);

% Rg 1.05 -----------------------------------------------------------------

Rg_105=zeros(5,1);
Rg_105_error=zeros(5,1);

N491051_Rg=importdata('N49_1.05_run1_Rg.txt');
N491052_Rg=importdata('N49_1.05_run2_Rg.txt');
N491053_Rg=importdata('N49_1.05_run3_Rg.txt');
N491054_Rg=importdata('N49_1.05_run4_Rg.txt');
N491055_Rg=importdata('N49_1.05_run5_Rg.txt');
N491051_Rg=N491051_Rg(L1:L2);
N491052_Rg=N491052_Rg(L1:L2);
N491053_Rg=N491053_Rg(L1:L2);
N491054_Rg=N491054_Rg(L1:L2);
N491055_Rg=N491055_Rg(L1:L2);

Rg_105(1)=mean([mean(N491051_Rg);mean(N491052_Rg);mean(N491053_Rg);mean(N491054_Rg);mean(N491055_Rg)]);
Rg_105_error(1)=std([mean(N491051_Rg);mean(N491052_Rg);mean(N491053_Rg);mean(N491054_Rg);mean(N491055_Rg)]);


N981051_Rg=importdata('N98_1.05_run1_Rg.txt');
N981052_Rg=importdata('N98_1.05_run2_Rg.txt');
N981053_Rg=importdata('N98_1.05_run3_Rg.txt');
N981054_Rg=importdata('N98_1.05_run4_Rg.txt');
N981055_Rg=importdata('N98_1.05_run5_Rg.txt');
N981051_Rg=N981051_Rg(L1:L2);
N981052_Rg=N981052_Rg(L1:L2);
N981053_Rg=N981053_Rg(L1:L2);
N981054_Rg=N981054_Rg(L1:L2);
N981055_Rg=N981055_Rg(L1:L2);

Rg_105(2)=mean([mean(N981051_Rg);mean(N981052_Rg);mean(N981053_Rg);mean(N981054_Rg);mean(N981055_Rg)]);
Rg_105_error(2)=std([mean(N981051_Rg);mean(N981052_Rg);mean(N981053_Rg);mean(N981054_Rg);mean(N981055_Rg)]);

NSP1051_Rg=importdata('NSP_1.05_run1_Rg.txt');
NSP1052_Rg=importdata('NSP_1.05_run2_Rg.txt');
NSP1053_Rg=importdata('NSP_1.05_run3_Rg.txt');
NSP1054_Rg=importdata('NSP_1.05_run4_Rg.txt');
NSP1055_Rg=importdata('NSP_1.05_run5_Rg.txt');
NSP1051_Rg=NSP1051_Rg(L1:L2);
NSP1052_Rg=NSP1052_Rg(L1:L2);
NSP1053_Rg=NSP1053_Rg(L1:L2);
NSP1054_Rg=NSP1054_Rg(L1:L2);
NSP1055_Rg=NSP1055_Rg(L1:L2);

Rg_105(3)=mean([mean(NSP1051_Rg);mean(NSP1052_Rg);mean(NSP1053_Rg);mean(NSP1054_Rg);mean(NSP1055_Rg)]);
Rg_105_error(3)=std([mean(NSP1051_Rg);mean(NSP1052_Rg);mean(NSP1053_Rg);mean(NSP1054_Rg);mean(NSP1055_Rg)]);

NUL1051_Rg=importdata('NUL_1.05_run1_Rg.txt');
NUL1052_Rg=importdata('NUL_1.05_run2_Rg.txt');
NUL1053_Rg=importdata('NUL_1.05_run3_Rg.txt');
NUL1054_Rg=importdata('NUL_1.05_run4_Rg.txt');
NUL1055_Rg=importdata('NUL_1.05_run5_Rg.txt');
NUL1051_Rg=NUL1051_Rg(L1:L2);
NUL1052_Rg=NUL1052_Rg(L1:L2);
NUL1053_Rg=NUL1053_Rg(L1:L2);
NUL1054_Rg=NUL1054_Rg(L1:L2);
NUL1055_Rg=NUL1055_Rg(L1:L2);

Rg_105(4)=mean([mean(NUL1051_Rg);mean(NUL1052_Rg);mean(NUL1053_Rg);mean(NUL1054_Rg);mean(NUL1055_Rg)]);
Rg_105_error(4)=std([mean(NUL1051_Rg);mean(NUL1052_Rg);mean(NUL1053_Rg);mean(NUL1054_Rg);mean(NUL1055_Rg)]);

NUS1051_Rg=importdata('NUS_1.05_run1_Rg.txt');
NUS1052_Rg=importdata('NUS_1.05_run2_Rg.txt');
NUS1053_Rg=importdata('NUS_1.05_run3_Rg.txt');
NUS1054_Rg=importdata('NUS_1.05_run4_Rg.txt');
NUS1055_Rg=importdata('NUS_1.05_run5_Rg.txt');
NUS1051_Rg=NUS1051_Rg(L1:L2);
NUS1052_Rg=NUS1052_Rg(L1:L2);
NUS1053_Rg=NUS1053_Rg(L1:L2);
NUS1054_Rg=NUS1054_Rg(L1:L2);
NUS1055_Rg=NUS1055_Rg(L1:L2);

Rg_105(5)=mean([mean(NUS1051_Rg);mean(NUS1052_Rg);mean(NUS1053_Rg);mean(NUS1054_Rg);mean(NUS1055_Rg)]);
Rg_105_error(5)=std([mean(NUS1051_Rg);mean(NUS1052_Rg);mean(NUS1053_Rg);mean(NUS1054_Rg);mean(NUS1055_Rg)]);


% FRET 1.05 -----------------------------------------------------------------

FRET_105=zeros(5,1);
FRET_105_error=zeros(5,1);

N491051_FRET=importdata('N49_1.05_run1_FRET.txt');
N491052_FRET=importdata('N49_1.05_run2_FRET.txt');
N491053_FRET=importdata('N49_1.05_run3_FRET.txt');
N491054_FRET=importdata('N49_1.05_run4_FRET.txt');
N491055_FRET=importdata('N49_1.05_run5_FRET.txt');
N491051_FRET=N491051_FRET(L1:L2);
N491052_FRET=N491052_FRET(L1:L2);
N491053_FRET=N491053_FRET(L1:L2);
N491054_FRET=N491054_FRET(L1:L2);
N491055_FRET=N491055_FRET(L1:L2);

FRET_105(1)=mean([mean(N491051_FRET);mean(N491052_FRET);mean(N491053_FRET);mean(N491054_FRET);mean(N491055_FRET)]);
FRET_105_error(1)=std([mean(N491051_FRET);mean(N491052_FRET);mean(N491053_FRET);mean(N491054_FRET);mean(N491055_FRET)]);


N981051_FRET=importdata('N98_1.05_run1_FRET.txt');
N981052_FRET=importdata('N98_1.05_run2_FRET.txt');
N981053_FRET=importdata('N98_1.05_run3_FRET.txt');
N981054_FRET=importdata('N98_1.05_run4_FRET.txt');
N981055_FRET=importdata('N98_1.05_run5_FRET.txt');
N981051_FRET=N981051_FRET(L1:L2);
N981052_FRET=N981052_FRET(L1:L2);
N981053_FRET=N981053_FRET(L1:L2);
N981054_FRET=N981054_FRET(L1:L2);
N981055_FRET=N981055_FRET(L1:L2);

FRET_105(2)=mean([mean(N981051_FRET);mean(N981052_FRET);mean(N981053_FRET);mean(N981054_FRET);mean(N981055_FRET)]);
FRET_105_error(2)=std([mean(N981051_FRET);mean(N981052_FRET);mean(N981053_FRET);mean(N981054_FRET);mean(N981055_FRET)]);

NSP1051_FRET=importdata('NSP_1.05_run1_FRET.txt');
NSP1052_FRET=importdata('NSP_1.05_run2_FRET.txt');
NSP1053_FRET=importdata('NSP_1.05_run3_FRET.txt');
NSP1054_FRET=importdata('NSP_1.05_run4_FRET.txt');
NSP1055_FRET=importdata('NSP_1.05_run5_FRET.txt');
NSP1051_FRET=NSP1051_FRET(L1:L2);
NSP1052_FRET=NSP1052_FRET(L1:L2);
NSP1053_FRET=NSP1053_FRET(L1:L2);
NSP1054_FRET=NSP1054_FRET(L1:L2);
NSP1055_FRET=NSP1055_FRET(L1:L2);

FRET_105(3)=mean([mean(NSP1051_FRET);mean(NSP1052_FRET);mean(NSP1053_FRET);mean(NSP1054_FRET);mean(NSP1055_FRET)]);
FRET_105_error(3)=std([mean(NSP1051_FRET);mean(NSP1052_FRET);mean(NSP1053_FRET);mean(NSP1054_FRET);mean(NSP1055_FRET)]);

NUL1051_FRET=importdata('NUL_1.05_run1_FRET.txt');
NUL1052_FRET=importdata('NUL_1.05_run2_FRET.txt');
NUL1053_FRET=importdata('NUL_1.05_run3_FRET.txt');
NUL1054_FRET=importdata('NUL_1.05_run4_FRET.txt');
NUL1055_FRET=importdata('NUL_1.05_run5_FRET.txt');
NUL1051_FRET=NUL1051_FRET(L1:L2);
NUL1052_FRET=NUL1052_FRET(L1:L2);
NUL1053_FRET=NUL1053_FRET(L1:L2);
NUL1054_FRET=NUL1054_FRET(L1:L2);
NUL1055_FRET=NUL1055_FRET(L1:L2);

FRET_105(4)=mean([mean(NUL1051_FRET);mean(NUL1052_FRET);mean(NUL1053_FRET);mean(NUL1054_FRET);mean(NUL1055_FRET)]);
FRET_105_error(4)=std([mean(NUL1051_FRET);mean(NUL1052_FRET);mean(NUL1053_FRET);mean(NUL1054_FRET);mean(NUL1055_FRET)]);

NUS1051_FRET=importdata('NUS_1.05_run1_FRET.txt');
NUS1052_FRET=importdata('NUS_1.05_run2_FRET.txt');
NUS1053_FRET=importdata('NUS_1.05_run3_FRET.txt');
NUS1054_FRET=importdata('NUS_1.05_run4_FRET.txt');
NUS1055_FRET=importdata('NUS_1.05_run5_FRET.txt');
NUS1051_FRET=NUS1051_FRET(L1:L2);
NUS1052_FRET=NUS1052_FRET(L1:L2);
NUS1053_FRET=NUS1053_FRET(L1:L2);
NUS1054_FRET=NUS1054_FRET(L1:L2);
NUS1055_FRET=NUS1055_FRET(L1:L2);

FRET_105(5)=mean([mean(NUS1051_FRET);mean(NUS1052_FRET);mean(NUS1053_FRET);mean(NUS1054_FRET);mean(NUS1055_FRET)]);
FRET_105_error(5)=std([mean(NUS1051_FRET);mean(NUS1052_FRET);mean(NUS1053_FRET);mean(NUS1054_FRET);mean(NUS1055_FRET)]);


% Rg 1.1 -----------------------------------------------------------------

Rg_11=zeros(5,1);
Rg_11_error=zeros(5,1);

N49111_Rg=importdata('N49_1.1_run1_Rg.txt');
N49112_Rg=importdata('N49_1.1_run2_Rg.txt');
N49113_Rg=importdata('N49_1.1_run3_Rg.txt');
N49114_Rg=importdata('N49_1.1_run4_Rg.txt');
N49115_Rg=importdata('N49_1.1_run5_Rg.txt');
N49111_Rg=N49111_Rg(L1:L2);
N49112_Rg=N49112_Rg(L1:L2);
N49113_Rg=N49113_Rg(L1:L2);
N49114_Rg=N49114_Rg(L1:L2);
N49115_Rg=N49115_Rg(L1:L2);

Rg_11(1)=mean([mean(N49111_Rg);mean(N49112_Rg);mean(N49113_Rg);mean(N49114_Rg);mean(N49115_Rg)]);
Rg_11_error(1)=std([mean(N49111_Rg);mean(N49112_Rg);mean(N49113_Rg);mean(N49114_Rg);mean(N49115_Rg)]);


N98111_Rg=importdata('N98_1.1_run1_Rg.txt');
N98112_Rg=importdata('N98_1.1_run2_Rg.txt');
N98113_Rg=importdata('N98_1.1_run3_Rg.txt');
N98114_Rg=importdata('N98_1.1_run4_Rg.txt');
N98115_Rg=importdata('N98_1.1_run5_Rg.txt');
N98111_Rg=N98111_Rg(L1:L2);
N98112_Rg=N98112_Rg(L1:L2);
N98113_Rg=N98113_Rg(L1:L2);
N98114_Rg=N98114_Rg(L1:L2);
N98115_Rg=N98115_Rg(L1:L2);

Rg_11(2)=mean([mean(N98111_Rg);mean(N98112_Rg);mean(N98113_Rg);mean(N98114_Rg);mean(N98115_Rg)]);
Rg_11_error(2)=std([mean(N98111_Rg);mean(N98112_Rg);mean(N98113_Rg);mean(N98114_Rg);mean(N98115_Rg)]);

NSP111_Rg=importdata('NSP_1.1_run1_Rg.txt');
NSP112_Rg=importdata('NSP_1.1_run2_Rg.txt');
NSP113_Rg=importdata('NSP_1.1_run3_Rg.txt');
NSP114_Rg=importdata('NSP_1.1_run4_Rg.txt');
NSP115_Rg=importdata('NSP_1.1_run5_Rg.txt');
NSP111_Rg=NSP111_Rg(L1:L2);
NSP112_Rg=NSP112_Rg(L1:L2);
NSP113_Rg=NSP113_Rg(L1:L2);
NSP114_Rg=NSP114_Rg(L1:L2);
NSP115_Rg=NSP115_Rg(L1:L2);

Rg_11(3)=mean([mean(NSP111_Rg);mean(NSP112_Rg);mean(NSP113_Rg);mean(NSP114_Rg);mean(NSP115_Rg)]);
Rg_11_error(3)=std([mean(NSP111_Rg);mean(NSP112_Rg);mean(NSP113_Rg);mean(NSP114_Rg);mean(NSP115_Rg)]);

NUL111_Rg=importdata('NUL_1.1_run1_Rg.txt');
NUL112_Rg=importdata('NUL_1.1_run2_Rg.txt');
NUL113_Rg=importdata('NUL_1.1_run3_Rg.txt');
NUL114_Rg=importdata('NUL_1.1_run4_Rg.txt');
NUL115_Rg=importdata('NUL_1.1_run5_Rg.txt');
NUL111_Rg=NUL111_Rg(L1:L2);
NUL112_Rg=NUL112_Rg(L1:L2);
NUL113_Rg=NUL113_Rg(L1:L2);
NUL114_Rg=NUL114_Rg(L1:L2);
NUL115_Rg=NUL115_Rg(L1:L2);

Rg_11(4)=mean([mean(NUL111_Rg);mean(NUL112_Rg);mean(NUL113_Rg);mean(NUL114_Rg);mean(NUL115_Rg)]);
Rg_11_error(4)=std([mean(NUL111_Rg);mean(NUL112_Rg);mean(NUL113_Rg);mean(NUL114_Rg);mean(NUL115_Rg)]);

NUS111_Rg=importdata('NUS_1.1_run1_Rg.txt');
NUS112_Rg=importdata('NUS_1.1_run2_Rg.txt');
NUS113_Rg=importdata('NUS_1.1_run3_Rg.txt');
NUS114_Rg=importdata('NUS_1.1_run4_Rg.txt');
NUS115_Rg=importdata('NUS_1.1_run5_Rg.txt');
NUS111_Rg=NUS111_Rg(L1:L2);
NUS112_Rg=NUS112_Rg(L1:L2);
NUS113_Rg=NUS113_Rg(L1:L2);
NUS114_Rg=NUS114_Rg(L1:L2);
NUS115_Rg=NUS115_Rg(L1:L2);

Rg_11(5)=mean([mean(NUS111_Rg);mean(NUS112_Rg);mean(NUS113_Rg);mean(NUS114_Rg);mean(NUS115_Rg)]);
Rg_11_error(5)=std([mean(NUS111_Rg);mean(NUS112_Rg);mean(NUS113_Rg);mean(NUS114_Rg);mean(NUS115_Rg)]);


% FRET 1.1 -----------------------------------------------------------------

FRET_11=zeros(5,1);
FRET_11_error=zeros(5,1);

N49111_FRET=importdata('N49_1.1_run1_FRET.txt');
N49112_FRET=importdata('N49_1.1_run2_FRET.txt');
N49113_FRET=importdata('N49_1.1_run3_FRET.txt');
N49114_FRET=importdata('N49_1.1_run4_FRET.txt');
N49115_FRET=importdata('N49_1.1_run5_FRET.txt');
N49111_FRET=N49111_FRET(L1:L2);
N49112_FRET=N49112_FRET(L1:L2);
N49113_FRET=N49113_FRET(L1:L2);
N49114_FRET=N49114_FRET(L1:L2);
N49115_FRET=N49115_FRET(L1:L2);

FRET_11(1)=mean([mean(N49111_FRET);mean(N49112_FRET);mean(N49113_FRET);mean(N49114_FRET);mean(N49115_FRET)]);
FRET_11_error(1)=std([mean(N49111_FRET);mean(N49112_FRET);mean(N49113_FRET);mean(N49114_FRET);mean(N49115_FRET)]);


N98111_FRET=importdata('N98_1.1_run1_FRET.txt');
N98112_FRET=importdata('N98_1.1_run2_FRET.txt');
N98113_FRET=importdata('N98_1.1_run3_FRET.txt');
N98114_FRET=importdata('N98_1.1_run4_FRET.txt');
N98115_FRET=importdata('N98_1.1_run5_FRET.txt');
N98111_FRET=N98111_FRET(L1:L2);
N98112_FRET=N98112_FRET(L1:L2);
N98113_FRET=N98113_FRET(L1:L2);
N98114_FRET=N98114_FRET(L1:L2);
N98115_FRET=N98115_FRET(L1:L2);

FRET_11(2)=mean([mean(N98111_FRET);mean(N98112_FRET);mean(N98113_FRET);mean(N98114_FRET);mean(N98115_FRET)]);
FRET_11_error(2)=std([mean(N98111_FRET);mean(N98112_FRET);mean(N98113_FRET);mean(N98114_FRET);mean(N98115_FRET)]);

NSP111_FRET=importdata('NSP_1.1_run1_FRET.txt');
NSP112_FRET=importdata('NSP_1.1_run2_FRET.txt');
NSP113_FRET=importdata('NSP_1.1_run3_FRET.txt');
NSP114_FRET=importdata('NSP_1.1_run4_FRET.txt');
NSP115_FRET=importdata('NSP_1.1_run5_FRET.txt');
NSP111_FRET=NSP111_FRET(L1:L2);
NSP112_FRET=NSP112_FRET(L1:L2);
NSP113_FRET=NSP113_FRET(L1:L2);
NSP114_FRET=NSP114_FRET(L1:L2);
NSP115_FRET=NSP115_FRET(L1:L2);

FRET_11(3)=mean([mean(NSP111_FRET);mean(NSP112_FRET);mean(NSP113_FRET);mean(NSP114_FRET);mean(NSP115_FRET)]);
FRET_11_error(3)=std([mean(NSP111_FRET);mean(NSP112_FRET);mean(NSP113_FRET);mean(NSP114_FRET);mean(NSP115_FRET)]);

NUL111_FRET=importdata('NUL_1.1_run1_FRET.txt');
NUL112_FRET=importdata('NUL_1.1_run2_FRET.txt');
NUL113_FRET=importdata('NUL_1.1_run3_FRET.txt');
NUL114_FRET=importdata('NUL_1.1_run4_FRET.txt');
NUL115_FRET=importdata('NUL_1.1_run5_FRET.txt');
NUL111_FRET=NUL111_FRET(L1:L2);
NUL112_FRET=NUL112_FRET(L1:L2);
NUL113_FRET=NUL113_FRET(L1:L2);
NUL114_FRET=NUL114_FRET(L1:L2);
NUL115_FRET=NUL115_FRET(L1:L2);

FRET_11(4)=mean([mean(NUL111_FRET);mean(NUL112_FRET);mean(NUL113_FRET);mean(NUL114_FRET);mean(NUL115_FRET)]);
FRET_11_error(4)=std([mean(NUL111_FRET);mean(NUL112_FRET);mean(NUL113_FRET);mean(NUL114_FRET);mean(NUL115_FRET)]);

NUS111_FRET=importdata('NUS_1.1_run1_FRET.txt');
NUS112_FRET=importdata('NUS_1.1_run2_FRET.txt');
NUS113_FRET=importdata('NUS_1.1_run3_FRET.txt');
NUS114_FRET=importdata('NUS_1.1_run4_FRET.txt');
NUS115_FRET=importdata('NUS_1.1_run5_FRET.txt');
NUS111_FRET=NUS111_FRET(L1:L2);
NUS112_FRET=NUS112_FRET(L1:L2);
NUS113_FRET=NUS113_FRET(L1:L2);
NUS114_FRET=NUS114_FRET(L1:L2);
NUS115_FRET=NUS115_FRET(L1:L2);

FRET_11(5)=mean([mean(NUS111_FRET);mean(NUS112_FRET);mean(NUS113_FRET);mean(NUS114_FRET);mean(NUS115_FRET)]);
FRET_11_error(5)=std([mean(NUS111_FRET);mean(NUS112_FRET);mean(NUS113_FRET);mean(NUS114_FRET);mean(NUS115_FRET)]);

% Plot data ----------------------------------------------------------------------------------------------------------------------

x=[1.0;20];
y=[1.0;20];

figure;
hold on;
plot(exp_FRET,FRET_09,'o','MarkerSize',24,'MarkerFaceColor',color1,'Linewidth',6,'Color',color1);
plot(exp_FRET,FRET_095,'o','MarkerSize',24,'MarkerFaceColor',color2,'Linewidth',6,'Color',color2);
plot(exp_FRET,FRET_10,'o','MarkerSize',24,'MarkerFaceColor',color4,'Linewidth',6,'Color',color4);
plot(exp_FRET,FRET_105,'o','MarkerSize',24,'MarkerFaceColor',color5,'Linewidth',6,'Color',color5);
plot(exp_FRET,FRET_11,'o','MarkerSize',24,'MarkerFaceColor',color6,'Linewidth',6,'Color',color6);
errorbar(exp_FRET,FRET_09,FRET_09_error,'o','MarkerSize',24,'MarkerFaceColor',color1,'Linewidth',6,'Color',color1);
errorbar(exp_FRET,FRET_095,FRET_095_error,'o','MarkerSize',24,'MarkerFaceColor',color2,'Linewidth',6,'Color',color2);
errorbar(exp_FRET,FRET_10,FRET_10_error,'o','MarkerSize',24,'MarkerFaceColor',color4,'Linewidth',6,'Color',color4);
errorbar(exp_FRET,FRET_105,FRET_105_error,'o','MarkerSize',24,'MarkerFaceColor',color5,'Linewidth',6,'Color',color5);
errorbar(exp_FRET,FRET_11,FRET_11_error,'o','MarkerSize',24,'MarkerFaceColor',color6,'Linewidth',6,'Color',color6);
errorbar(exp_FRET,FRET_09,exp_FRET_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color1,'Linewidth',6,'Color',color1);
errorbar(exp_FRET,FRET_095,exp_FRET_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color2,'Linewidth',6,'Color',color2);
errorbar(exp_FRET,FRET_10,exp_FRET_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color4,'Linewidth',6,'Color',color4);
errorbar(exp_FRET,FRET_105,exp_FRET_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color5,'Linewidth',6,'Color',color5);
errorbar(exp_FRET,FRET_11,exp_FRET_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color6,'Linewidth',6,'Color',color6);
plot(x,y,'--k','MarkerSize',24,'Linewidth',6);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4);
legend({'0.9','0.95','1.0','1.05','1.1'},'location','southeast');
axis([2 11 2 11]);
box on;

figure;
hold on;
plot(exp_Rg,Rg_09,'o','MarkerSize',24,'MarkerFaceColor',color1,'Linewidth',6,'Color',color1);
plot(exp_Rg,Rg_095,'o','MarkerSize',24,'MarkerFaceColor',color2,'Linewidth',6,'Color',color2);
plot(exp_Rg,Rg_10,'o','MarkerSize',24,'MarkerFaceColor',color4,'Linewidth',6,'Color',color4);
plot(exp_Rg,Rg_105,'o','MarkerSize',24,'MarkerFaceColor',color5,'Linewidth',6,'Color',color5);
plot(exp_Rg,Rg_11,'o','MarkerSize',24,'MarkerFaceColor',color6,'Linewidth',6,'Color',color6);
errorbar(exp_Rg,Rg_09,Rg_09_error,'o','MarkerSize',24,'MarkerFaceColor',color1,'Linewidth',6,'Color',color1);
errorbar(exp_Rg,Rg_095,Rg_095_error,'o','MarkerSize',24,'MarkerFaceColor',color2,'Linewidth',6,'Color',color2);
errorbar(exp_Rg,Rg_10,Rg_10_error,'o','MarkerSize',24,'MarkerFaceColor',color4,'Linewidth',6,'Color',color4);
errorbar(exp_Rg,Rg_105,Rg_105_error,'o','MarkerSize',24,'MarkerFaceColor',color5,'Linewidth',6,'Color',color5);
errorbar(exp_Rg,Rg_11,Rg_11_error,'o','MarkerSize',24,'MarkerFaceColor',color6,'Linewidth',6,'Color',color6);
errorbar(exp_Rg,Rg_09,exp_Rg_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color1,'Linewidth',6,'Color',color1);
errorbar(exp_Rg,Rg_095,exp_Rg_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color2,'Linewidth',6,'Color',color2);
errorbar(exp_Rg,Rg_10,exp_Rg_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color4,'Linewidth',6,'Color',color4);
errorbar(exp_Rg,Rg_105,exp_Rg_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color5,'Linewidth',6,'Color',color5);
errorbar(exp_Rg,Rg_11,exp_Rg_error,'o',"horizontal",'MarkerSize',24,'MarkerFaceColor',color6,'Linewidth',6,'Color',color6);
plot(x,y,'--k','MarkerSize',24,'Linewidth',6);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4);
legend({'0.9','0.95','1.0','1.05','1.1'},'location','southeast');
axis([1.0 5 1.0 5]);
box on;

% Calculate chi^2 ----------------------------------------------------------------------------------------------------------------------

chi_Rg_09 = (exp_Rg-Rg_09).^2 ./ (exp_Rg_error).^2;
chi_FRET_09 = (exp_FRET-FRET_09).^2 ./ (exp_FRET_error).^2;
chi_Rg_095 = (exp_Rg-Rg_095).^2 ./ (exp_Rg_error).^2;
chi_FRET_095 = (exp_FRET-FRET_095).^2 ./ (exp_FRET_error).^2;
chi_Rg_10 = (exp_Rg-Rg_10).^2 ./ (exp_Rg_error).^2;
chi_FRET_10 = (exp_FRET-FRET_10).^2 ./ (exp_FRET_error).^2;
chi_Rg_105 = (exp_Rg-Rg_105).^2 ./ (exp_Rg_error).^2;
chi_FRET_105 = (exp_FRET-FRET_105).^2 ./ (exp_FRET_error).^2;
chi_Rg_11 = (exp_Rg-Rg_11).^2 ./ (exp_Rg_error).^2;
chi_FRET_11 = (exp_FRET-FRET_11).^2 ./ (exp_FRET_error).^2;

chi_Rg_mean=[mean(chi_Rg_09);mean(chi_Rg_095);mean(chi_Rg_10);mean(chi_Rg_105);mean(chi_Rg_11)];

chi_FRET_mean=[mean(chi_FRET_09);mean(chi_FRET_095);mean(chi_FRET_10);mean(chi_FRET_105);mean(chi_FRET_11)];

chi_combine=mean([chi_FRET_mean,chi_Rg_mean],2);

x=[0.9:0.05:1.1];

figure;
hold on;
plot(x,chi_combine,'-ok','MarkerSize',24,'Linewidth',6,'MarkerFaceColor',[0,0,0]);
set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',4,'yscale','log');
%set(gca,'FontSize',52,'FontName','Helvetica','Linewidth',6);
axis([0.9 1.1 8 25])
box on;
