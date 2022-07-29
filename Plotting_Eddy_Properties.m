clear; clc; close all;
% CMEMS eddy tracking output
cmems_rep = '/Volumes/NEMO3/GOFFISH/Emily/Programs/Eddy_Tracking_CMEMS_NEW/EXTRACTION/';
% Anticyclonic Eddies
load([cmems_rep 'EDDY_PROPERTIES/AE_Properties.mat']);
Amp_AE = 100*double(Amplitude); % amplitude in cm
EKE_AE = 10000*double(EKE);    % EKE in cm2 s-2
Radius_AE = double(Radius);
Xcenter=double(Xcenter); Ycenter=double(Ycenter); 
e2=EKE_AE; 
e2(~isnan(EKE_AE))=1; 
Nanti=sum(e2,1,'omitnan'); 
% Cyclonic Eddies
load([cmems_rep 'EDDY_PROPERTIES/CE_Properties.mat']);
Amp_CE = 100*double(Amplitude);
EKE_CE = 10000*double(EKE);
Radius_CE = double(Radius);
date_num2=double(date_num);datelabels=datestr(date_num2);
Xcenter=double(Xcenter); Ycenter=double(Ycenter); 
e2=EKE_CE; 
e2(~isnan(EKE_CE))=1; 
Ncyclo=sum(e2,1,'omitnan'); 
% Average all eddies per day
Radius_CE=median(Radius_CE,'omitnan'); Radius_AE=median(Radius_AE,'omitnan'); 
Amp_CE=median(Amp_CE,'omitnan'); Amp_AE=median(Amp_AE,'omitnan'); 
EKE_CE=median(EKE_CE,'omitnan'); EKE_AE=median(EKE_AE,'omitnan'); 
Num_AE=Nanti; Num_CE=Ncyclo;
% Rename variables and extract 2019-2021
Amp_AE_CMEMS = Amp_AE(8036:9131);Amp_CE_CMEMS = Amp_CE(8036:9131);
EKE_AE_CMEMS = EKE_AE(8036:9131);EKE_CE_CMEMS = EKE_CE(8036:9131);
Num_AE_CMEMS = Num_AE(8036:9131);Num_CE_CMEMS = Num_CE(8036:9131);
Radius_AE_CMEMS = Radius_AE(8036:9131);Radius_CE_CMEMS = Radius_CE(8036:9131); 
% Calculate moving average to smooth data
Amp_AE_CMEMS=movmean(Amp_AE_CMEMS,30,'omitnan');
Amp_CE_CMEMS=movmean(Amp_CE_CMEMS,30,'omitnan');
Num_AE_CMEMS=movmean(Num_AE_CMEMS,30,'omitnan'); 
Num_CE_CMEMS=movmean(Num_CE_CMEMS,30,'omitnan'); 
Radius_AE_CMEMS=movmean(Radius_AE_CMEMS,30,'omitnan');
Radius_CE_CMEMS=movmean(Radius_CE_CMEMS,30,'omitnan');
EKE_AE_CMEMS=movmean(EKE_AE_CMEMS,30,'omitnan');
EKE_CE_CMEMS=movmean(EKE_CE_CMEMS,30,'omitnan');
%% HYCOM eddy tracking output
% HYCOM eddy tracking output
hycom_rep = '/Volumes/Kate-Research/Data/Eddy_Extraction/';
% Anticyclonic Eddies
load([hycom_rep 'EDDY_PROPERTIES/AE_Properties.mat']);
Amp_AE = 100*double(Amplitude); % amplitude in cm
EKE_AE = 10000*double(EKE);    % EKE in cm2 s-2
Radius_AE = double(Radius);
Xcenter=double(Xcenter); Ycenter=double(Ycenter); 
e2=EKE_AE; 
e2(~isnan(EKE_AE))=1; 
Nanti=sum(e2,1,'omitnan'); 
% Cyclonic Eddies
load([hycom_rep 'EDDY_PROPERTIES/CE_Properties.mat']);
Amp_CE = 100*double(Amplitude);
EKE_CE = 10000*double(EKE);
Radius_CE = double(Radius);
date_num2=double(date_num);datelabels=datestr(date_num2);
Xcenter=double(Xcenter); Ycenter=double(Ycenter); 
e2=EKE_CE; 
e2(~isnan(EKE_CE))=1; 
Ncyclo=sum(e2,1,'omitnan'); 
% Average all eddies per day
Radius_CE=median(Radius_CE,'omitnan'); Radius_AE=median(Radius_AE,'omitnan'); 
Amp_CE=median(Amp_CE,'omitnan'); Amp_AE=median(Amp_AE,'omitnan'); 
EKE_CE=median(EKE_CE,'omitnan'); EKE_AE=median(EKE_AE,'omitnan'); 
Num_AE=Nanti; Num_CE=Ncyclo;
% Rename variables
Amp_AE_HYCOM = Amp_AE; Amp_CE_HYCOM = Amp_CE;
EKE_AE_HYCOM = EKE_AE; EKE_CE_HYCOM = EKE_CE;
Num_AE_HYCOM = Num_AE; Num_CE_HYCOM = Num_CE;
Radius_AE_HYCOM = Radius_AE; Radius_CE_HYCOM = Radius_CE; 
% Calculate moving average to smooth data
Amp_AE_HYCOM=movmean(Amp_AE_HYCOM,30,'omitnan');
Amp_CE_HYCOM=movmean(Amp_CE_HYCOM,30,'omitnan');
Num_AE_HYCOM=movmean(Num_AE_HYCOM,30,'omitnan'); 
Num_CE_HYCOM=movmean(Num_CE_HYCOM,30,'omitnan'); 
Radius_AE_HYCOM=movmean(Radius_AE_HYCOM,30,'omitnan');
Radius_CE_HYCOM=movmean(Radius_CE_HYCOM,30,'omitnan');
EKE_AE_HYCOM=movmean(EKE_AE_HYCOM,30,'omitnan');
EKE_CE_HYCOM=movmean(EKE_CE_HYCOM,30,'omitnan');
%% Plotting Code
addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/Corinne_toolbox');
addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/m_map');
addpath('/Volumes/Lacie-SAN/SAN2/MATLAB/mytoolbox');% Sarahs toolbox :P
addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/cmocean_v2.0/cmocean');

% cd /Volumes/NEMO3/GOFFISH/Emily/Data/EddyPropertiesMatFiles
% load Eddy_PropertiesHYCOM.mat
% Amp_AE_HYCOM = Amp_AE;Amp_CE_HYCOM = Amp_CE;
% EKE_AE_HYCOM = EKE_AE;EKE_CE_HYCOM = EKE_CE;
% Num_AE_HYCOM = Num_AE;Num_CE_HYCOM = Num_CE;
% Radius_AE_HYCOM = Radius_AE;Radius_CE_HYCOM = Radius_CE; 

% cd /Volumes/NEMO3/GOFFISH/Emily/Data/EddyPropertiesMatFiles
% load Eddy_PropertiesCMEMS.mat
% Amp_AE_CMEMS = Amp_AE(8036:9131);Amp_CE_CMEMS = Amp_CE(8036:9131);
% EKE_AE_CMEMS = EKE_AE(8036:9131);EKE_CE_CMEMS = EKE_CE(8036:9131);
% Num_AE_CMEMS = Num_AE(8036:9131);Num_CE_CMEMS = Num_CE(8036:9131);
% Radius_AE_CMEMS = Radius_AE(8036:9131);Radius_CE_CMEMS = Radius_CE(8036:9131); 

cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load Ocolor_eddies_avg.mat
CMEMS_AE = movmean(bigmeanCMEMS_AE,30,'omitnan');
CMEMS_CE = movmean(bigmeanCMEMS_CE,30,'omitnan');
HYCOM_AE = movmean(bigmeanHYCOM_AE,30,'omitnan');
HYCOM_CE = movmean(bigmeanHYCOM_CE,30,'omitnan');

X = linspace(1, 1096, 1096);
fs = 22;
figure('Position',[1 1 2300 1300])
ha=tight_subplot(5,1,[.02 .015],[.065 .04],.095);

axes(ha(1))
plot(X, Amp_AE_HYCOM,'color',[255/256, 153/256, 51/256],'linewidth',2.5);
hold on
plot(X, Amp_AE_CMEMS,'color','r','linewidth',2.5);
hold on
plot(X, Amp_CE_HYCOM,'color',[0, 204/256, 255/256],'linewidth',2.5);
hold on
plot(X, Amp_CE_CMEMS,'color','b','linewidth',2.5);
set(gca,'XTick', 1:30.4444:1096)
set(gca,'XTickLabels', {''})
xlim([1 1096])
ylim([2 15]);
ylabel('Eddy Amplitude (cm)')
l1 = legend('AE HYCOM','AE Altimetry','CE HYCOM','CE Altimetry');
legend('Location','northwest')
legend('Orientation','Horizontal')
grid on
set(gca,'fontsize',fs)

axes(ha(2))
plot(X, Num_AE_HYCOM,'color',[255/256, 153/256, 51/256],'linewidth',2.5);
hold on
plot(X, Num_AE_CMEMS,'color','r','linewidth',2.5);
hold on
plot(X, Num_CE_HYCOM,'color',[0, 204/256, 255/256],'linewidth',2.5);
hold on
plot(X, Num_CE_CMEMS,'color','b','linewidth',2.5);
set(gca,'XTick', 1:30.4444:1096)
set(gca,'XTickLabels', {''})
xlim([1 1096])
ylim([2 20])
set(gca,'YTick', 0:5:20)
ylabel('Number of Eddies')
grid on
set(gca,'fontsize',fs)


axes(ha(3))
plot(X, Radius_AE_HYCOM,'color',[255/256, 153/256, 51/256],'linewidth',2.5);
hold on
plot(X, Radius_AE_CMEMS,'color','r','linewidth',2.5);
hold on
plot(X, Radius_CE_HYCOM,'color',[0, 204/256, 255/256],'linewidth',2.5);
hold on
plot(X, Radius_CE_CMEMS,'color','b','linewidth',2.5);
set(gca,'XTick', 1:30.4444:1096)
set(gca,'XTickLabels', {''})
xlim([1 1096])
ylim([50 120])
ylabel('Eddy Radius (km)')
grid on
set(gca,'fontsize',fs)

axes(ha(4))
plot(X, EKE_AE_HYCOM,'color',[255/256, 153/256, 51/256],'linewidth',2.5);
hold on
plot(X, EKE_AE_CMEMS,'color','r','linewidth',2.5);
hold on
plot(X, EKE_CE_HYCOM,'color',[0, 204/256, 255/256],'linewidth',2.5);
hold on
plot(X, EKE_CE_CMEMS,'color','b','linewidth',2.5);
set(gca,'XTick', 1:30.4444:1096)
set(gca,'YTick', 0:200:800)
set(gca,'XTickLabels', {''})
xlim([1 1096])
ylim([0 800]);
ylabel('Eddy EKE (cm^2 s^-^2)')
grid on
set(gca,'fontsize',fs)

axes(ha(5))
plot(X, HYCOM_AE,'color',[255/256, 153/256, 51/256],'linewidth',2.5);
hold on
plot(X, CMEMS_AE,'color','r','linewidth',2.5);
hold on
plot(X, HYCOM_CE,'color',[0, 204/256, 255/256],'linewidth',2.5);
hold on
plot(X, CMEMS_CE,'color','b','linewidth',2.5);
set(gca,'XTick', 1:30.4444:1096)
set(gca,'YTick', -.04:.02:.04);
xlim([1 1096])
ylim([-.04 .04]);
ylabel('Chl-a (mg/m^3)')
% set(gca, 'YScale', 'log')
grid on
set(gca,'fontsize',fs)
% set(gca,'XTickLabels', {'2019','2020','2021',''});
% hold on;
% set(gca,'XTick', 1:30.5:365)
set(gca,'XTickLabels', {'2019','F','M','A','M','J','J','A','S','O','N','D','2020','F','M','A','M','J','J','A','S','O','N','D','2021','F','M','A','M','J','J','A','S','O','N','D'})
xtickangle(45)
% 
% axes(ha(5))
% plot(X, EKE_LC97,'color',[255/256, 140/256, 0],'linewidth',2.5);
% hold on
% plot(X, EKE_LC98,'color',[255/256, 215/256, 0],'linewidth',2.5);
% hold on
% plot(X, EKE_LC14,'color',[67/256, 110/256, 238/256],'linewidth',2.5);
% hold on
% plot(X, EKE_LC15,'color',[85/256,26/256,139/256],'linewidth',2.5);
% xlim([1 365])
% ylabel('LC EKE (cm^2 s^-^2)')
% grid on
% set(gca,'fontsize',fs,'XTickLabelRotation',-45)


cd('/Volumes/Kate-Research/Figures')
saveas(gca,'Tseries_5Panel_19to21_CMEMS_HYCOM.tif');