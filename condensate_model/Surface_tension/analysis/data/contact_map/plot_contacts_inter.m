clear all;
close all;

mymap=[43   140 190
          78   179 211
          123  204 196
          168  221 181
          204  235 197
          224  243 219
          255  255 255
          254  232 200
          253  212 158
          253  187 132
          252  141 89
          239  101 72
          230  80  60]/255;

mymap2=zeros(100,3);
mymap2=mymap2+1;



for i=1:100
    if i < 51
        index=i-1;
        mymap2(i,1)=mymap2(i,1)-index*0.0015;
        mymap2(i,2)=mymap2(i,2)-index*0.0015;
    else
        index=i-50;
        mymap2(i,1)=mymap2(i,1)-(index*0.0185+50*0.0015);
        mymap2(i,2)=mymap2(i,2)-(index*0.0185+50*0.0015);
    end
end

% Plot Nup62c -----------------------------------------------------------

Nup62c_inter1=importdata('Nup62_Nup58_Nup54_32_inter.txt');
Nup62c_inter2=importdata('Nup62_Nup58_Nup54_32_2_inter.txt');
Nup62c_inter3=importdata('Nup62_Nup58_Nup54_32_3_inter.txt');
Nup62c_inter4=importdata('Nup62_Nup58_Nup54_32_4_inter.txt');
Nup62c_inter5=importdata('Nup62_Nup58_Nup54_32_5_inter.txt');
Nup62c_inter6=importdata('Nup62_Nup58_Nup54_32_6_inter.txt');
Nup62c_inter7=importdata('Nup62_Nup58_Nup54_32_7_inter.txt');
Nup62c_inter8=importdata('Nup62_Nup58_Nup54_32_8_inter.txt');
Nup62c_inter9=importdata('Nup62_Nup58_Nup54_32_9_inter.txt');
Nup62c_inter10=importdata('Nup62_Nup58_Nup54_32_10_inter.txt');
Nup62c_inter11=importdata('Nup62_Nup58_Nup54_32_11_inter.txt');
Nup62c_inter12=importdata('Nup62_Nup58_Nup54_32_12_inter.txt');
Nup62c_inter13=importdata('Nup62_Nup58_Nup54_32_13_inter.txt');
Nup62c_inter14=importdata('Nup62_Nup58_Nup54_32_14_inter.txt');
Nup62c_inter15=importdata('Nup62_Nup58_Nup54_32_15_inter.txt');

% only simulations with one condensate
Nup62c_inter_v1=(Nup62c_inter1+Nup62c_inter3+Nup62c_inter4+Nup62c_inter5+Nup62c_inter6+Nup62c_inter7+Nup62c_inter9+Nup62c_inter10+Nup62c_inter11+Nup62c_inter12)./10;


N=max(size(Nup62c_inter_v1));

Nup62c_inter_v2=zeros(N,N);

for i=1:N
    Nup62c_inter_v2(i,i)=2*Nup62c_inter_v1(i,i);
    for j=i+1:N
        Nup62c_inter_v2(j,i)=(Nup62c_inter_v1(i,j)+Nup62c_inter_v1(j,i));
        Nup62c_inter_v2(i,j)=Nup62c_inter_v2(j,i);
    end
end

fig=figure('Renderer', 'painters', 'Position', [0 0 750 750]);
hm = heatmap(Nup62c_inter_v2);
colormap(flipud(parula));
set(gca,'FontSize',36,'FontName','Helvetica');
caxis([0 0.01]);

 
% Get underlying axis handle
origState = warning('query', 'MATLAB:structOnObject');
cleanup = onCleanup(@()warning(origState));
warning('off','MATLAB:structOnObject')
S = struct(hm); % Undocumented
ax = S.Axes;    % Undocumented
clear('cleanup')
% Remove grids
hm.GridVisible = 'off';
% Place lines around selected columns and row
% Assumes columns and rows are 1 unit in size!
col = [1, 111, 508, 753, 933, 1107, 1436, 1629];    
row = [1, 111, 508, 753, 933, 1107, 1436, 1629];
xline(ax, [col-.5], 'k-','LineWidth',1.5); % see footnotes [1,2]
yline(ax, [row-.5], 'k-','LineWidth',1.5); % see footnotes [1,2]
print(fig,'inter_contacts.png','-dpng');
