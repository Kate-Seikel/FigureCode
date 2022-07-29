%% Matlab &&
% Kate Seikel 07/15/22
clc; clear variables;
% addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/Corinne_toolbox');
% addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/m_map'); % Don't use this
addpath('/Volumes/Lacie-SAN/SAN2/MATLAB/mytoolbox');% Sarahs toolbox :P
% addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/cmocean_v2.0/cmocean');
addpath('/Volumes/Kate-Research/MatlabPrograms');
%% CMEMS DATA
main_rep = '/Volumes/Kate-Research/MatlabPrograms/Eddy_Tracking_CMEMS/EXTRACTION/EDDY_PROPERTIES';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANTICYCLONINC CMEMS DATA
% Ocolor data
load('/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/Ocolor_eddies_avg_spatial.mat');
% AE Properties
load([main_rep '/AE_Properties.mat'])
Amp_AE = 100*double(Amplitude); % amplitude in cm
EKE_AE = 10000*double(EKE);    % EKE in cm2 s-2
Radius_AE = double(Radius);
Xcenter_AE = double(Xcenter);Ycenter_AE = double(Ycenter); Xcenter_AE = Xcenter_AE+360;
Xcentroid_AE = double(Xcentroid);Ycentroid_AE = double(Ycentroid); Xcentroid_AE = Xcentroid_AE+360;
% SST and SSS Data(TBD)
load('/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SST_eddies_avg.mat');
load('/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SSS_eddies_avg.mat');
% % RENAME to _sat - not separating by temp anomaly
Amp_AE_sat=Amp_AE;
SSTA_AE_sat=SSTbigmeanCMEMS_AE;
Rad_AE_sat=Radius_AE;
NV_AE_sat=bigmeanCMEMS_AE;
Xcenter_AE_sat=Xcenter_AE;
Xcentroid_AE_sat=Xcentroid_AE;
Ycenter_AE_sat=Ycenter_AE;
Ycentroid_AE_sat=Ycentroid_AE;
SSSA_AE_sat=SSSbigmeanCMEMS_AE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CYCLONIC CMEMS DATA
% CE Properties
load([main_rep '/CE_Properties.mat'])
Amp_CE = 100*double(Amplitude);
EKE_CE = 10000*double(EKE);
Radius_CE = double(Radius);
Xcenter_CE = double(Xcenter);Ycenter_CE = double(Ycenter); Xcenter_CE = Xcenter_CE+360;
Xcentroid_CE = double(Xcentroid);Ycentroid_CE = double(Ycentroid); Xcentroid_CE = Xcentroid_CE+360;
% SST and SSS Data(TBD)
% RENAME to _sat - not separating by temp anomaly
Amp_CE_sat=Amp_CE;
SSTA_CE_sat=SSTbigmeanCMEMS_CE;
Rad_CE_sat=Radius_CE;
NV_CE_sat=bigmeanCMEMS_CE;
Xcenter_CE_sat=Xcenter_CE;
Xcentroid_CE_sat=Xcentroid_CE;
Ycenter_CE_sat=Ycenter_CE;
Ycentroid_CE_sat=Ycentroid_CE;
SSSA_CE_sat=SSSbigmeanCMEMS_CE;
%%
Xc = [Xcentroid_AE(:);Xcentroid_CE(:)];
Yc = [Ycentroid_AE(:);Ycentroid_CE(:)];
%% CALCUL DES MOYENNES DES PROPRIÉTÉS (Rayon, Amplitude, EKE) SUR UNE GRILLE SPATIALE
Resol = 1;   % resolution de la grille en degrés.
Xg = min([Xcentroid_AE_sat(:);Xcentroid_CE_sat(:)])-Resol/2:Resol:max([Xcentroid_AE_sat(:);Xcentroid_CE_sat(:)])+Resol/2;
Yg = min([Ycentroid_AE_sat(:);Ycentroid_CE_sat(:)])-Resol/2:Resol:max([Ycentroid_AE_sat(:);Ycentroid_CE_sat(:)])+Resol/2;

ind = find(isfinite(Xcentroid_AE_sat)==1);
N_AElayergrid = ac_SumOnGrid(Xcentroid_AE_sat(ind(3)),Ycentroid_AE_sat(ind(3)),ones(1,length(ind(1))),Xg(:),Yg(:)); % finds the grid location of the first eddy
for i=ind(2:end)
    %ind=i;
    try
        N_AElayer = ac_SumOnGrid(Xcentroid_AE_sat(i),Ycentroid_AE_sat(i),ones(1,length(i)),Xg(:),Yg(:)); % finds the grid location of the first eddy
          N_AElayergrid=N_AElayer+N_AElayergrid;

    catch
        continue
    end
end
N_AE_sat=N_AElayergrid;
R_AE_sat = ac_SumOnGrid(Xcentroid_AE_sat(ind),Ycentroid_AE_sat(ind),Rad_AE_sat(ind),Xg(:),Yg(:))./N_AE_sat; R_AE_sat(R_AE_sat>10000)=NaN;
A_AE_sat = ac_SumOnGrid(Xcentroid_AE_sat(ind),Ycentroid_AE_sat(ind),Amp_AE_sat(ind),Xg(:),Yg(:))./N_AE_sat; A_AE_sat(A_AE_sat>10000)=NaN;
NV_AE_sat(isnan(NV_AE_sat))=0;
NV_AE_sat = ac_SumOnGrid(Xcentroid_AE_sat(ind),Ycentroid_AE_sat(ind),NV_AE_sat(ind),Xg(:),Yg(:))./N_AE_sat; % NV_AE_sat(NV_AE_sat>100)=NaN;
SSTA_AE_sat(isnan(SSTA_AE_sat))=0;
SSTA_AE_sat = ac_SumOnGrid(Xcentroid_AE_sat(ind),Ycentroid_AE_sat(ind),SSTA_AE_sat(ind),Xg(:),Yg(:))./N_AE_sat; SSTA_AE_sat(SSTA_AE_sat>100)=NaN;
SSSA_AE_sat(isnan(SSSA_AE_sat))=0;
SSSA_AE_sat = ac_SumOnGrid(Xcentroid_AE_sat(ind),Ycentroid_AE_sat(ind),SSSA_AE_sat(ind),Xg(:),Yg(:))./N_AE_sat; % SSSA_AE_sat(SSSA_AE_sat>100)=NaN;

%Resol = 2;   % resolution de la grille en degrés.
Xg2 = min([Xcentroid_CE_sat(:);Xcentroid_CE_sat(:)])-Resol/2:Resol:max([Xcentroid_CE_sat(:);Xcentroid_CE_sat(:)])+Resol/2;
Yg2 = min([Ycentroid_CE_sat(:);Ycentroid_CE_sat(:)])-Resol/2:Resol:max([Ycentroid_CE_sat(:);Ycentroid_CE_sat(:)])+Resol/2;

ind = find(isfinite(Xcentroid_CE_sat)==1);
N_CElayergrid = ac_SumOnGrid(Xcentroid_CE_sat(ind(1)),Ycentroid_CE_sat(ind(1)),ones(1,length(ind(1))),Xg2(:),Yg2(:)); % finds the grid location of the first eddy
for i=ind(2:end)
    %ind=i;
    try
        N_CElayer = ac_SumOnGrid(Xcentroid_CE_sat(i),Ycentroid_CE_sat(i),ones(1,length(i)),Xg2(:),Yg2(:)); % finds the grid location of the first eddy
          N_CElayergrid=N_CElayer+N_CElayergrid;

    catch
        continue
    end
end
N_CE_sat=N_CElayergrid;
R_CE_sat = ac_SumOnGrid(Xcentroid_CE_sat(ind),Ycentroid_CE_sat(ind),Rad_CE_sat(ind),Xg2(:),Yg2(:))./N_CE_sat; R_CE_sat(R_CE_sat>10000)=NaN;
A_CE_sat = ac_SumOnGrid(Xcentroid_CE_sat(ind),Ycentroid_CE_sat(ind),Amp_CE_sat(ind),Xg2(:),Yg2(:))./N_CE_sat; A_CE_sat(A_CE_sat>10000)=NaN;
NV_CE_sat(isnan(NV_CE_sat))=0;
NV_CE_sat = ac_SumOnGrid(Xcentroid_CE_sat(ind),Ycentroid_CE_sat(ind),NV_CE_sat(ind),Xg2(:),Yg2(:))./N_CE_sat; % NV_CE_sat(NV_CE_sat>100)=NaN;
SSTA_CE_sat(isnan(SSTA_CE_sat))=0;
SSTA_CE_sat = ac_SumOnGrid(Xcentroid_CE_sat(ind),Ycentroid_CE_sat(ind),SSTA_CE_sat(ind),Xg2(:),Yg2(:))./N_CE_sat; SSTA_CE_sat(SSTA_CE_sat>100)=NaN;
SSSA_CE_sat(isnan(SSSA_CE_sat))=0;
SSSA_CE_sat = ac_SumOnGrid(Xcentroid_CE_sat(ind),Ycentroid_CE_sat(ind),SSSA_CE_sat(ind),Xg2(:),Yg2(:))./N_CE_sat; SSSA_CE_sat(SSSA_CE_sat>100)=NaN;

%% HYCOM (mod) data
%%
main_rep = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES';
load([main_rep '/AE_Properties.mat'])
Amp_AE = 100*double(Amplitude); % amplitude in cm
EKE_AE = 10000*double(EKE);    % EKE in cm2 s-2
Radius_AE = double(Radius);
Xcenter_AE = double(Xcenter);Ycenter_AE = double(Ycenter);  Xcenter_AE = Xcenter_AE+360;
Xcentroid_AE = double(Xcentroid);Ycentroid_AE = double(Ycentroid); Xcentroid_AE = Xcentroid_AE+360;
% load([main_rep '/OcolorEddiesHYCOM.mat']);
load('/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SST_eddies_avg.mat');
load('/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SSS_eddies_avg.mat');
% RENAME to _mod - not separating by temp anomaly
Amp_AE_mod=Amp_AE;
SSTA_AE_mod=SSTbigmeanHYCOM_AE;
Rad_AE_mod=Radius_AE;
NV_AE_mod=bigmeanHYCOM_AE;
Xcenter_AE_mod=Xcenter_AE;
Xcentroid_AE_mod=Xcentroid_AE;
Ycenter_AE_mod=Ycenter_AE;
Ycentroid_AE_mod=Ycentroid_AE;
SSSA_AE_mod=SSSbigmeanHYCOM_AE;

%
load([main_rep '/CE_Properties.mat'])
Amp_CE = 100*double(Amplitude);
EKE_CE = 10000*double(EKE);
Radius_CE = double(Radius);
Xcenter_CE = double(Xcenter);Ycenter_CE = double(Ycenter); Xcenter_CE = Xcenter_CE+360;
Xcentroid_CE = double(Xcentroid);Ycentroid_CE = double(Ycentroid); Xcentroid_CE = Xcentroid_CE+360;
% load([main_rep '/CE_Properties_SSTA_SSSA_dec5_ssh']);

% RENAME to _mod - not separating by temp anomaly
Amp_CE_mod=Amp_CE;
SSTA_CE_mod=SSTbigmeanHYCOM_CE;
Rad_CE_mod=Radius_CE;
NV_CE_mod=bigmeanHYCOM_CE;
Xcenter_CE_mod=Xcenter_CE;
Xcentroid_CE_mod=Xcentroid_CE;
Ycenter_CE_mod=Ycenter_CE;
Ycentroid_CE_mod=Ycentroid_CE;
SSSA_CE_mod=SSSbigmeanHYCOM_CE;

%%
Xc = [Xcentroid_AE(:);Xcentroid_CE(:)];
Yc = [Ycentroid_AE(:);Ycentroid_CE(:)];

%% CALCUL DES MOYENNES DES PROPRIÉTÉS (Rayon, Amplitude, EKE) SUR UNE GRILLE SPATIALE
Resol = 1;   % resolution de la grille en degrés.
Xg3 = min([Xcentroid_AE_mod(:);Xcentroid_CE_mod(:)])-Resol/2:Resol:max([Xcentroid_AE_mod(:);Xcentroid_CE_mod(:)])+Resol/2;
Yg3 = min([Ycentroid_AE_mod(:);Ycentroid_CE_mod(:)])-Resol/2:Resol:max([Ycentroid_AE_mod(:);Ycentroid_CE_mod(:)])+Resol/2;

ind = find(isfinite(Xcentroid_AE_mod)==1);
N_AElayergrid = ac_SumOnGrid(Xcentroid_AE_mod(ind(3)),Ycentroid_AE_mod(ind(3)),ones(1,length(ind(1))),Xg3(:),Yg3(:)); % finds the grid location of the first eddy
for i=ind(2:end)
    %ind=i;
    try
        N_AElayer = ac_SumOnGrid(Xcentroid_AE_mod(i),Ycentroid_AE_mod(i),ones(1,length(i)),Xg3(:),Yg3(:)); % finds the grid location of the first eddy
          N_AElayergrid=N_AElayer+N_AElayergrid;

    catch
        continue
    end
end
N_AE_mod=N_AElayergrid;
R_AE_mod = ac_SumOnGrid(Xcentroid_AE_mod(ind),Ycentroid_AE_mod(ind),Rad_AE_mod(ind),Xg3(:),Yg3(:))./N_AE_mod; R_AE_mod(R_AE_mod>10000)=NaN;
A_AE_mod = ac_SumOnGrid(Xcentroid_AE_mod(ind),Ycentroid_AE_mod(ind),Amp_AE_mod(ind),Xg3(:),Yg3(:))./N_AE_mod; A_AE_mod(A_AE_mod>10000)=NaN;
NV_AE_mod(isnan(NV_AE_mod)) = 0;
NV_AE_mod = ac_SumOnGrid(Xcentroid_AE_mod(ind),Ycentroid_AE_mod(ind),NV_AE_mod(ind),Xg3(:),Yg3(:))./N_AE_mod; % NV_AE_mod=NV_AE_mod./10^8; NV_AE_mod(NV_AE_mod>10000000)=NaN;
SSTA_AE_mod(isnan(SSTA_AE_mod)) = 0;
SSTA_AE_mod = ac_SumOnGrid(Xcentroid_AE_mod(ind),Ycentroid_AE_mod(ind),SSTA_AE_mod(ind),Xg3(:),Yg3(:))./N_AE_mod;  SSTA_AE_mod(SSTA_AE_mod>100)=NaN;
SSSA_AE_mod(isnan(SSSA_AE_mod)) = 0;
SSSA_AE_mod = ac_SumOnGrid(Xcentroid_AE_mod(ind),Ycentroid_AE_mod(ind),SSSA_AE_mod(ind),Xg3(:),Yg3(:))./N_AE_mod; % SSSA_AE_mod(SSSA_AE_mod>100)=NaN;

%Resol = 2;   % resolution de la grille en degrés.
Xg4 = min([Xcentroid_CE_mod(:);Xcentroid_CE_mod(:)])-Resol/2:Resol:max([Xcentroid_CE_mod(:);Xcentroid_CE_mod(:)])+Resol/2;
Yg4 = min([Ycentroid_CE_mod(:);Ycentroid_CE_mod(:)])-Resol/2:Resol:max([Ycentroid_CE_mod(:);Ycentroid_CE_mod(:)])+Resol/2;

ind = find(isfinite(Xcentroid_CE_mod)==1);
N_CElayergrid = ac_SumOnGrid(Xcentroid_CE_mod(ind(1)),Ycentroid_CE_mod(ind(1)),ones(1,length(ind(1))),Xg4(:),Yg4(:)); % finds the grid location of the first eddy
for i=ind(2:end)
    %ind=i;
    try
        N_CElayer = ac_SumOnGrid(Xcentroid_CE_mod(i),Ycentroid_CE_mod(i),ones(1,length(i)),Xg4(:),Yg4(:)); % finds the grid location of the first eddy
          N_CElayergrid=N_CElayer+N_CElayergrid;

    catch
        continue
    end
end
N_CE_mod=N_CElayergrid;
R_CE_mod = ac_SumOnGrid(Xcentroid_CE_mod(ind),Ycentroid_CE_mod(ind),Rad_CE_mod(ind),Xg4(:),Yg4(:))./N_CE_mod; R_CE_mod(R_CE_mod>10000)=NaN;
A_CE_mod = ac_SumOnGrid(Xcentroid_CE_mod(ind),Ycentroid_CE_mod(ind),Amp_CE_mod(ind),Xg4(:),Yg4(:))./N_CE_mod; A_CE_mod(A_CE_mod>10000)=NaN;
NV_CE_mod(isnan(NV_CE_mod)) = 0;
NV_CE_mod = ac_SumOnGrid(Xcentroid_CE_mod(ind),Ycentroid_CE_mod(ind),NV_CE_mod(ind),Xg4(:),Yg4(:))./N_CE_mod; % NV_CE_mod=NV_CE_mod./10^8; NV_CE_mod(NV_CE_mod>10000000)=NaN;
SSTA_CE_mod(isnan(SSTA_CE_mod)) = 0;
SSTA_CE_mod = ac_SumOnGrid(Xcentroid_CE_mod(ind),Ycentroid_CE_mod(ind),SSTA_CE_mod(ind),Xg4(:),Yg4(:))./N_CE_mod;  SSTA_CE_mod(SSTA_CE_mod>100)=NaN;
SSSA_CE_mod(isnan(SSSA_CE_mod)) = 0;
SSSA_CE_mod = ac_SumOnGrid(Xcentroid_CE_mod(ind),Ycentroid_CE_mod(ind),SSSA_CE_mod(ind),Xg4(:),Yg4(:))./N_CE_mod;  SSSA_CE_mod(SSSA_CE_mod>100)=NaN;
%% FIGURES

% Number of Eddies
%-------------------


m_proj('mercator','lon',[260 280],'lat',[18 32])
figure('Position',[1 1 1600 900])

ha=tight_subplot(4,6,[.01 .01],[.12 .06],[.15 .07]);
%subplot(2,1,1)
N_AE_sat(N_AE_sat<=5)=NaN;
axes(ha(1))
my_m_pcolor(Xg,Yg,N_AE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(a)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([0 200])
set(gca,'fontsize',23)
title('Num of Eddies')
%ylabel('Satellite')
m_text(252,22,'Satellite','fontsize', 25, 'rotation',90)

N_AE_mod(N_AE_mod<=5)=NaN;
axes(ha(7))
my_m_pcolor(Xg3,Yg3,N_AE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[]);
caxis([0 200])
set(gca,'fontsize',23)%,'Position',[0.415 0.7502 0.3825 0.2348])
colormap jet
m_text(260.5,30.5,'(g)','color','w','fontsize',35)
m_text(252,22,'Model','fontsize', 25, 'rotation',90)
m_text(246,27,'Anticyclonic','fontsize', 30, 'rotation',90)

N_CE_sat(N_CE_sat<=5)=NaN;
axes(ha(13))
%N_CE_sat(isnan(N_CE_sat))=0;
my_m_pcolor(Xg2,Yg2,N_CE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(m)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([0 200])
set(gca,'fontsize',23)
%ylabel('Satellite')
m_text(252,22,'Satellite','fontsize', 25, 'rotation',90)

xlab1='100^o W';

N_CE_mod(N_CE_mod<=5)=NaN;
axes(ha(19))
my_m_pcolor(Xg4,Yg4,N_CE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xtick',[264 272 280]);
caxis([0 200])
colormap jet
set(gca,'fontsize',23)%,'Position',[0.415 0.7502 0.3825 0.2348])
 %ylabel('Model')
m_text(252,22,'Model','fontsize', 25, 'rotation',90)
m_text(246,27,'Cyclonic','fontsize', 30, 'rotation',90)

m_text(260.5,30.5,'(s)','color','w','fontsize',35)
 col1=colorbar;
 set(col1,'location','southoutside')
 set(col1,'Position',[ 0.155    0.05    0.1    0.02])

% Radius 

R_AE_sat(isnan(N_AE_sat))=NaN;
axes(ha(2))
my_m_pcolor(Xg,Yg,R_AE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(b)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([40 160])
set(gca,'fontsize',23)
title('Radius (km)')

R_AE_mod(isnan(N_AE_mod))=NaN;
axes(ha(8))
my_m_pcolor(Xg3,Yg3,R_AE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(h)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([40 160])
set(gca,'fontsize',23)

R_CE_sat(isnan(N_CE_sat))=NaN;
axes(ha(14))
my_m_pcolor(Xg2,Yg2,R_CE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(n)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([40 160])
set(gca,'fontsize',23)

R_CE_mod(isnan(N_CE_mod))=NaN;
axes(ha(20))
my_m_pcolor(Xg4,Yg4,R_CE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(t)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([40 160])
set(gca,'fontsize',23)
 col2=colorbar;
 set(col2,'location','southoutside')
 set(col2,'Position',[ 0.29    0.05    0.1    0.02])


% Amplitude
A_AE_sat(isnan(N_AE_sat))=NaN;
axes(ha(3))
my_m_pcolor(Xg,Yg,A_AE_sat')
hold on
shading flat
hold on
m_line([48 63], [4 4],'color','w','linewidth',4)
m_line([63 63], [4 25],'color','w','linewidth',4)

m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(c)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([0 30])
set(gca,'fontsize',23)
title('Amplitude (cm)')

A_AE_mod(isnan(N_AE_mod))=NaN;
axes(ha(9))
my_m_pcolor(Xg3,Yg3,A_AE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(i)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([0 30])
set(gca,'fontsize',23)

A_CE_sat(isnan(N_CE_sat))=NaN;
axes(ha(15))
my_m_pcolor(Xg2,Yg2,A_CE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(o)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([0 30])
set(gca,'fontsize',23)

A_CE_mod(isnan(N_CE_mod))=NaN;
axes(ha(21))
my_m_pcolor(Xg4,Yg4,A_CE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(u)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([0 30])
set(gca,'fontsize',23)
 col3=colorbar;
 set(col3,'location','southoutside')
 set(col3,'Position',[ 0.42    0.05    0.1    0.02])

% Rossby Numb / Chl-a
% CHLA_AE_sat(isnan(N_AE_sat))=NaN;
NV_AE_sat(isnan(N_AE_sat))=NaN;
axes(ha(4))
% my_m_pcolor(Xg(2:19),Yg,CHLA_AE_sat')
my_m_pcolor(Xg,Yg,NV_AE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(d)','color','w','fontsize',35)
% c1=colorbar; set(c1,'location','northoutside','position',[0.5450 0.72 0.1217 0.02]);
%title('Number of CEs')
caxis([-.04 .04])
set(gca,'fontsize',23)
title('Chl-a') % (x10^{-3})')

NV_AE_mod(isnan(N_AE_mod))=NaN;
axes(ha(10))
my_m_pcolor(Xg3,Yg3,NV_AE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(j)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.04 .04])
set(gca,'fontsize',23)

% CHLA_CE_sat(isnan(N_CE_sat))=NaN;
NV_CE_sat(isnan(N_CE_sat))=NaN;
axes(ha(16))
% my_m_pcolor(Xg(2:19),Yg,CHLA_CE_sat')
my_m_pcolor(Xg2,Yg2,NV_CE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(p)','color','w','fontsize',35)
%title('Number of CEs')
caxis([-.04 .04])
set(gca,'fontsize',23)

% delete(ha(22))
NV_CE_mod(isnan(N_CE_mod))=NaN;
axes(ha(22))
my_m_pcolor(Xg4,Yg4,NV_CE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(v)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.04 .04])
set(gca,'fontsize',23)
 col4=colorbar;
 set(col4,'location','southoutside')
 set(col4,'Position',[ 0.555    0.05    0.1    0.02])

%% SSTA 
SSTA_AE_sat(isnan(N_AE_sat))=NaN;
axes(ha(5))
my_m_pcolor(Xg,Yg,SSTA_AE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(e)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.85 .85])
set(gca,'fontsize',23)
title('SSTA (^oC)')

SSTA_AE_mod(isnan(N_AE_mod))=NaN;
axes(ha(11))
my_m_pcolor(Xg3,Yg3,SSTA_AE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(k)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.85 .85])
set(gca,'fontsize',23)

SSTA_CE_sat(isnan(N_CE_sat))=NaN;
axes(ha(17))
my_m_pcolor(Xg2,Yg2,SSTA_CE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(q)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.85 .85])
set(gca,'fontsize',23)

SSTA_CE_mod(isnan(N_CE_mod))=NaN;
axes(ha(23))
my_m_pcolor(Xg4,Yg4,SSTA_CE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(w)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.85 .85])
set(gca,'fontsize',23)
 col5=colorbar;
 set(col5,'location','southoutside')
 set(col5,'Position',[ 0.69    0.05    0.1    0.02])
 
 % SSSA 
SSSA_AE_sat(isnan(N_AE_sat))=NaN;
axes(ha(6))
my_m_pcolor(Xg,Yg,SSSA_AE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(f)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.4 .41])
set(gca,'fontsize',23)
title('SSSA (psu)')
% m_text(282,31,'Surface Int.','fontsize', 25, 'rotation',270)

SSSA_AE_mod(isnan(N_AE_mod))=NaN;
axes(ha(12))
my_m_pcolor(Xg3,Yg3,SSSA_AE_mod')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(l)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.4 .41])
set(gca,'fontsize',23)
%m_text(282,33,'Subsurface Int.','fontsize', 25, 'rotation',270)

SSSA_CE_sat(isnan(N_CE_sat))=NaN;
axes(ha(18))
my_m_pcolor(Xg2,Yg2,SSSA_CE_sat')
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'xticklabel',[],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(r)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.4 .41])
set(gca,'fontsize',23)
%m_text(282,33,'Subsurface Int.','fontsize', 25, 'rotation',270)

SSSA_CE_mod(isnan(N_CE_mod))=NaN;
axes(ha(24))
my_m_pcolor(Xg4,Yg4,SSSA_CE_mod');
hold on
shading flat
hold on
m_coast('patch',[.6 .6 .6]);
m_grid('box','fancy','fontsize',18,'xtick',[264 272 280],'yticklabel',[])%,'YAxisLocation','right');
m_text(260.5,30.5,'(x)','color','w','fontsize',35)
%colorbar
%title('Number of CEs')
caxis([-.4 .41])
set(gca,'fontsize',23)
%m_text(282,31,'Surface Int.','fontsize', 25, 'rotation',270)
 col6=colorbar;
 set(col6,'location','southoutside')
 set(col6,'Position',[ 0.825    0.05    0.1    0.02])