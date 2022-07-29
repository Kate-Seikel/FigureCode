%% Matlab &&
% Kate's code :P 06/29/22 
% This code is to create files to use in in the
% spatial map. I'm gonna pull the AE/CE contours and use those to make a
% file of average SSTA and SSSA(A is for anomaly) per eddy. 
% This is only for panels e,f,k,l,q,r,w,x,d,j,p,v
% Use already created CMEMS and HYCOM Eddy data for a,b,c,g,h,i,m,n,o,s,t,u
clear variables; clc;
addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/Corinne_toolbox');
addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/m_map'); % Don't use this
addpath('/Volumes/Lacie-SAN/SAN2/MATLAB/mytoolbox');% Sarahs toolbox :P
addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/cmocean_v2.0/cmocean');
addpath('/Volumes/Kate-Research/MatlabPrograms');
%% Paths to Data
% SATELLITE DATA!!!!!!!!! (SSS & SST)
SSTdatapath = '/Volumes/Lacie-SAN/SAN2/CMC_SST/CMC-L4-GLOB-v3.0-1deg'; % Path to SST Data
years = ["2019",'2020',"2021"]; % List of years
%% PART ONE Satellite Data
bigHYSST = NaN(141,181,1096); % 1096 is the total # of days from 19-21 
filenumSST = 1;
for i = 1:length(years)
    SSTpath = (SSTdatapath);
    yrstr=num2str(years(i));
    pathtoyearly = [SSTpath '/' yrstr '/'];
    cd(pathtoyearly);
    list = dir('*.nc');
    filenamesSST = string({list.name});
    filenamesSST = filenamesSST';
    daysinfileSST = length(filenamesSST);
    for j = 1:daysinfileSST
        % Read in adt data
        lat = ncread(filenamesSST(j),'lat'); 
        lon = ncread(filenamesSST(j),'lon'); 
        sst = ncread(filenamesSST(j),'analysed_sst'); 
        lon = double(lon); lat = double(lat); 
        lamin=nearestpoint(18, lat); lamax=nearestpoint(32,lat);
        lomin=nearestpoint(-98,lon); lomax=nearestpoint(-80,lon);
        plotlonSST = lon(lomin:lomax); plotlatSST = lat(lamin:lamax);
        [plotlonSST,plotlatSST] = meshgrid(plotlonSST,plotlatSST);
        % Narrow array down to the Gulf of Mexico
        plotSST = sst(lomin:lomax,lamin:lamax);
        plotSST = plotSST';
        % Create array to hold data to average later
        bigHYSST(:,:,filenumSST) = plotSST;
        filenumSST = filenumSST + 1;
    end
end
%% Saving SST file
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/SST_sat.mat';
save(filename_out,'bigHYSST','plotlatSST','plotlonSST');
clear variables; clc;
cd /Volumes/Kate-Research/Data/Eddy_Extraction/Variables
load('SST_sat.mat');
%% Compute anomalies
SSTsatref = NaN(141,181,365); % Array to hold the climatological mean
for i = 1:366 % Leap day is in 2020
    meanholder = NaN(141,181,3);
    if i==60 % LEAP DAY!! Calculate avg of Feb 28/March1st of 2019/2021 to get climate mean for the leap day
        lpdayholder = NaN(141,181,2);
        lpdayholder(:,:,1) = bigHYSST(:,:,59);
        lpdayholder(:,:,2) = bigHYSST(:,:,60);
        mean19 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,1) = mean19;
        meanholder(:,:,2) = bigHYSST(:,:,425);
        lpdayholder = NaN(141,181,2);
        lpdayholder(:,:,1) = bigHYSST(:,:,790);
        lpdayholder(:,:,2) = bigHYSST(:,:,791);
        mean21 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,3) = mean21;
        lpdaymeanSST = mean(meanholder,3,'omitnan');
    elseif i<60
        meanholder(:,:,1) = bigHYSST(:,:,i);
        meanholder(:,:,2) = bigHYSST(:,:,(i+365));
        meanholder(:,:,3) = bigHYSST(:,:,(i+731));
        SSTsatref(:,:,i) = mean(meanholder,3,'omitnan');
    elseif i>60
        meanholder(:,:,1) = bigHYSST(:,:,i-1);
        meanholder(:,:,2) = bigHYSST(:,:,(i+365));
        meanholder(:,:,3) = bigHYSST(:,:,(i+730));
        SSTsatref(:,:,i-1) = mean(meanholder,3,'omitnan');
    end
end
SSTanom = NaN(141,181,1096);
for i = 1:1096
    if i == 425 % Leap Day!
        SSTanom(:,:,i) = bigHYSST(:,:,i)-lpdaymeanSST;
    elseif i<=365
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i);
    elseif i>=365 && i<425
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i-365);
    elseif i>425 && i<=731
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i-366);
    elseif i>=732 && i<=1096
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i-731);
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/SSTanom.mat';
save(filename_out,'SSTanom');
%% HYCOM Eddies
cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load('Eddy_Contours.mat');
SSTanom = SSTanom;
AE_cont = AEs; CE_cont = CEs;
l = length(AE_cont(:,1,1));
SST_eddys_AE = cell(l,1096);
% Pull out SST for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    AE_cont_x = AE_cont(:,i,1); AE_cont_y = AE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        xa = double(AE_cont_x{j}); ya = double(AE_cont_y{j});
        SST = SSTanom(:,:,i);
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonSST,plotlatSST,xa,ya);% Take only datapoints within eddy boundary
%             outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
%             outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST and subtract the mean of the outside points to
            % get rid of seasonal variations
            SST_eddys_AE(j,i) = {SST};
        else 
            disp(j)
            disp(i)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(CE_cont(:,1,1));
SST_eddys_CE = cell(l,1096);
% Pull out SST for cyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    CE_cont_x = CE_cont(:,i,1); CE_cont_y = CE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        SST = SSTanom(:,:,i);
        xa = double(CE_cont_x{j}); ya = double(CE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonSST,plotlatSST,xa,ya);% Take only datapoints within eddy boundary
%             outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
%             outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST 
            SST_eddys_CE(j,i) = {SST};
        else 
            disp(j)
            disp(i)
        end
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SSTEddiesHYCOM.mat';
save(filename_out,'SST_eddys_AE','SST_eddys_CE','plotlonSST','plotlatSST','-v7.3');
% CMEMS Eddies
cd /Volumes/Kate-Research/MatlabPrograms/Eddy_Tracking_CMEMS/EXTRACTION/EDDY_PROPERTIES
load('Eddy_Contours.mat');
cd /Volumes/Kate-Research/Data/Eddy_Extraction
AE_cont = AEs; CE_cont = CEs;
l = length(AE_cont(:,1,1));
SST_eddys_AE = cell(l,1096);
% Pull out SST for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    AE_cont_x = AE_cont(:,i,1); AE_cont_y = AE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        SST = SSTanom(:,:,i);
        xa = double(AE_cont_x{j}); ya = double(AE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonSST,plotlatSST,xa,ya);% Take only datapoints within eddy boundary
            outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
            outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST 
            SST_eddys_AE(j,i) = {SST};
        else 
            disp(j)
            disp(i)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(CE_cont(:,1,1));
SST_eddys_CE = cell(l,1096);
% Pull out SST for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    CE_cont_x = CE_cont(:,i,1); CE_cont_y = CE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        SST = SSTanom(:,:,i);
        xa = double(CE_cont_x{j}); ya = double(CE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonSST,plotlatSST,xa,ya);% Take only datapoints within eddy boundary
            outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
            outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST 
            SST_eddys_CE(j,i) = {SST};
        else 
            disp(j)
            disp(i)
        end
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SSTEddiesCMEMS.mat';
save(filename_out,'SST_eddys_AE','SST_eddys_CE','plotlonSST','plotlatSST','-v7.3');
clear; clc;
%% Take time average of CMEMS/HYCOM for CEs and AEs
clear; clc;
cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load('SSTEddiesHYCOM.mat');
% HYCOM AEs
mean1 = NaN(21,1);
SSTbigmeanHYCOM_AE = NaN(21,1096);
for i = 1:1096
    AE = SST_eddys_AE(:,i);
    for j = 1:13
        if isempty(SST_eddys_AE(j,i))==0
            % take average SST from said eddy
            AE1 = double(AE{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1(j) = AE_mean;
        else
        end
    end
%     meannew = median(mean1,'omitnan');
    SSTbigmeanHYCOM_AE(:,i) = mean1;
end
% HYCOM CEs
mean1 = NaN(21,1);
SSTbigmeanHYCOM_CE = NaN(21,1096);
for i = 1:1096
    CE = SST_eddys_CE(:,i);
    for j = 1:18
        if isempty(SST_eddys_CE(j,i))==0
            % take average SST from said eddy
            CE1 = double(CE{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1(j) = CE_mean;
        else
        end
    end
%     meannew = median(mean1,'omitnan');
    SSTbigmeanHYCOM_CE(:,i) = mean1;
end
% Do the same for CMEMS
clearvars -except SSTbigmeanHYCOM_CE SSTbigmeanHYCOM_AE
load('SSTEddiesCMEMS.mat');
% HYCOM AEs
mean1 = NaN(21,1);
SSTbigmeanCMEMS_AE = NaN(21,1096);
for i = 1:1096
    AE = SST_eddys_AE(:,i);
    for j = 1:13
        if isempty(SST_eddys_AE(j,i))==0
            % take average SST from said eddy
            AE1 = double(AE{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1(j) = AE_mean;
        else
        end
    end
%     meannew = median(mean1,'omitnan');
    SSTbigmeanCMEMS_AE(:,i) = mean1;
end
% HYCOM CEs
mean1 = NaN(21,1);
SSTbigmeanCMEMS_CE = NaN(21,1096);
for i = 1:1096
    CE = SST_eddys_CE(:,i);
    for j = 1:18
        if isempty(SST_eddys_CE(j,i))==0
            % take average SST from said eddy
            CE1 = double(CE{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1(j) = CE_mean;
        else
        end
    end
%     meannew = median(mean1,'omitnan');
    SSTbigmeanCMEMS_CE(:,i) = mean1;
end
% climanom_CMEMS_CE = mean(bigmeanCMEMS_CE,'all','omitnan');
% climanom_CMEMS_AE = mean(bigmeanCMEMS_AE,'all','omitnan');
% climanom_HYCOM_CE = mean(bigmeanHYCOM_CE,'all','omitnan');
% climanom_HYCOM_AE = mean(bigmeanHYCOM_AE,'all','omitnan');
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SST_eddies_avg.mat';
save(filename_out,'SSTbigmeanCMEMS_CE','SSTbigmeanCMEMS_AE','SSTbigmeanHYCOM_AE','SSTbigmeanHYCOM_CE','-v7.3');
%% Do the same for satellite Salinity
clear; clc;
pathto2019 = '/Volumes/Lacie-SAN/SAN2/SMOS_BEC_GLOBAL/v2_L4/daily/2019';
pathtoSSSdata = '/Volumes/Lacie-SAN/SAN2/SMAP-RSS/SMAP-SSSV5.0/L3/8day_running'; % 2010-2021
years = ["2019",'2020',"2021"]; % List of years
filenumSSS = 1;
bigSSS = NaN(57,73,1096);
for i = 1:length(years) % cycles through the years!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yearnum = str2double(years(i));
    Altpath = (pathtoSSSdata);
    yrstr=num2str(years(i));
    pathtoyearly = [Altpath '/' yrstr '/'];
    cd(pathtoyearly);
    list = dir('*.nc');
    filenamesSSS2 = string({list.name});
    filenamesSSS2 = filenamesSSS2';
    if yearnum == 2019
        cd(pathto2019);
        list = dir('*.nc');
        filenamesSSS = string({list.name});
        filenamesSSS = filenamesSSS';
        daysinfileSSS = length(filenamesSSS);
        for j=1:365
            if j>=168 && j<= 206
                cd(pathto2019);
                lat = ncread(filenamesSSS(j),'lat');
                lon = ncread(filenamesSSS(j),'lon');
                sss = ncread(filenamesSSS(j),'sss');
                lon = double(lon); lat = double(lat); 
                lamin=nearestpoint(18, lat); lamax=nearestpoint(32,lat);
                lomin=nearestpoint(-98,lon); lomax=nearestpoint(-80,lon);
                plotlonSSSmiss = lon(lomin:lomax); plotlatSSSmiss = lat(lamin:lamax);
                [plotlonSSSmiss,plotlatSSSmiss] = meshgrid(plotlonSSSmiss,plotlatSSSmiss);
                plotSSS = sss(lomin:lomax,lamin:lamax);
                plotSSS = plotSSS';
                % Interpolate to match the other SSS files that aren't
                % missing
                plotlonSSSmissint = interp2(plotlonSSSmiss,plotlatSSSmiss,plotlonSSSmiss,plotlonSSS,plotlatSSS);
                plotlatSSSmissint = interp2(plotlonSSSmiss,plotlatSSSmiss,plotlatSSSmiss,plotlonSSS,plotlatSSS);
                plotSSS = interp2(plotlonSSSmiss,plotlatSSSmiss,plotSSS,plotlonSSS,plotlatSSS);
                bigSSS(:,:,filenumSSS) = plotSSS;
                disp(filenumSSS);
                filenumSSS = filenumSSS + 1;
            elseif j>=1 && j<168
                cd(pathtoyearly);
                lat = ncread(filenamesSSS2(j),'lat');
                lon = ncread(filenamesSSS2(j),'lon');
                sss = ncread(filenamesSSS2(j),'sss_smap');
                lon = double(lon); lat = double(lat); 
                % the lon is on 0 - 360 instead of -180 - 180
                dat1=lon(1:720);
                dat2=lon(721:end);
                dat3=dat2-360;
                dat4=NaN(1440,1);
                dat4(1:720)=dat3;
                dat4(721:end)=dat1;
                lon=dat4;
                % Rearrange Data to fit new longitudes
                u1=sss(1:720,:);
                u2=sss(721:end,:);
                u3=NaN(1440,720);
                u3(1:720,:)=u2;
                u3(721:end,:)=u1;
                % Narrow array down to the Gulf of Mexico
                lamin=nearestpoint(18, lat); lamax=nearestpoint(32,lat);
                lomin=nearestpoint(-98,lon); lomax=nearestpoint(-80,lon);
                plotlonSSS = lon(lomin:lomax); plotlatSSS = lat(lamin:lamax);
                [plotlonSSS,plotlatSSS] = meshgrid(plotlonSSS,plotlatSSS);
                plotSSS = u3(lomin:lomax,lamin:lamax);
                plotSSS = plotSSS';
                bigSSS(:,:,filenumSSS) = plotSSS;
                disp(filenumSSS);
                filenumSSS = filenumSSS + 1;
            elseif j > 206
                cd(pathtoyearly);
                lat = ncread(filenamesSSS2(j-39),'lat');
                lon = ncread(filenamesSSS2(j-39),'lon');
                sss = ncread(filenamesSSS2(j-39),'sss_smap');
                lon = double(lon); lat = double(lat); 
                % the lon is on 0 - 360 instead of -180 - 180
                dat1=lon(1:720);
                dat2=lon(721:end);
                dat3=dat2-360;
                dat4=NaN(1440,1);
                dat4(1:720)=dat3;
                dat4(721:end)=dat1;
                lon=dat4;
                % Rearrange Data to fit new longitudes
                u1=sss(1:720,:);
                u2=sss(721:end,:);
                u3=NaN(1440,720);
                u3(1:720,:)=u2;
                u3(721:end,:)=u1;
                % Narrow array down to the Gulf of Mexico
                lamin=nearestpoint(18, lat); lamax=nearestpoint(32,lat);
                lomin=nearestpoint(-98,lon); lomax=nearestpoint(-80,lon);
                plotlonSSS = lon(lomin:lomax); plotlatSSS = lat(lamin:lamax);
                [plotlonSSS,plotlatSSS] = meshgrid(plotlonSSS,plotlatSSS);
                plotSSS = u3(lomin:lomax,lamin:lamax);
                plotSSS = plotSSS';
                bigSSS(:,:,filenumSSS) = plotSSS;
                disp(filenumSSS);
                filenumSSS = filenumSSS + 1;
            end
        end
    elseif yearnum == 2020
        for j = 1:366
            if j==155
                holder = NaN(1440,720,2);
                lat = ncread(filenamesSSS2(j),'lat');
                lon = ncread(filenamesSSS2(j),'lon');
                holder(:,:,1) = ncread(filenamesSSS2(j-1),'sss_smap');
                holder(:,:,2) = ncread(filenamesSSS2(j+1),'sss_smap');
                sss155 = mean(holder,3,'omitnan');
                lon = double(lon); lat = double(lat); 
                % the lon is on 0 - 360 instead of -180 - 180
                dat1=lon(1:720);
                dat2=lon(721:end);
                dat3=dat2-360;
                dat4=NaN(1440,1);
                dat4(1:720)=dat3;
                dat4(721:end)=dat1;
                lon=dat4;
                % Rearrange Data to fit new longitudes
                u1=sss155(1:720,:);
                u2=sss155(721:end,:);
                u3=NaN(1440,720);
                u3(1:720,:)=u2;
                u3(721:end,:)=u1;
                % Narrow array down to the Gulf of Mexico
                lamin=nearestpoint(18, lat); lamax=nearestpoint(32,lat);
                lomin=nearestpoint(-98,lon); lomax=nearestpoint(-80,lon);
                plotlonSSS = lon(lomin:lomax); plotlatSSS = lat(lamin:lamax);
                [plotlonSSS,plotlatSSS] = meshgrid(plotlonSSS,plotlatSSS);
                plotSSS = u3(lomin:lomax,lamin:lamax);
                plotSSS = plotSSS';
                bigSSS(:,:,filenumSSS) = plotSSS;
                disp(filenumSSS);
                filenumSSS = filenumSSS + 1;
            elseif j<155
                lat = ncread(filenamesSSS2(j),'lat');
                lon = ncread(filenamesSSS2(j),'lon');
                sss = ncread(filenamesSSS2(j),'sss_smap');
                lon = double(lon); lat = double(lat); 
                % the lon is on 0 - 360 instead of -180 - 180
                dat1=lon(1:720);
                dat2=lon(721:end);
                dat3=dat2-360;
                dat4=NaN(1440,1);
                dat4(1:720)=dat3;
                dat4(721:end)=dat1;
                lon=dat4;
                % Rearrange Data to fit new longitudes
                u1=sss(1:720,:);
                u2=sss(721:end,:);
                u3=NaN(1440,720);
                u3(1:720,:)=u2;
                u3(721:end,:)=u1;
                % Narrow array down to the Gulf of Mexico
                lamin=nearestpoint(18, lat); lamax=nearestpoint(32,lat);
                lomin=nearestpoint(-98,lon); lomax=nearestpoint(-80,lon);
                plotlonSSS = lon(lomin:lomax); plotlatSSS = lat(lamin:lamax);
                [plotlonSSS,plotlatSSS] = meshgrid(plotlonSSS,plotlatSSS);
                plotSSS = u3(lomin:lomax,lamin:lamax);
                plotSSS = plotSSS';
                bigSSS(:,:,filenumSSS) = plotSSS;
                disp(filenumSSS);
                filenumSSS = filenumSSS + 1;
            elseif j>155
                lat = ncread(filenamesSSS2(j-1),'lat');
                lon = ncread(filenamesSSS2(j-1),'lon');
                sss = ncread(filenamesSSS2(j-1),'sss_smap');
                lon = double(lon); lat = double(lat); 
                % the lon is on 0 - 360 instead of -180 - 180
                dat1=lon(1:720);
                dat2=lon(721:end);
                dat3=dat2-360;
                dat4=NaN(1440,1);
                dat4(1:720)=dat3;
                dat4(721:end)=dat1;
                lon=dat4;
                % Rearrange Data to fit new longitudes
                u1=sss(1:720,:);
                u2=sss(721:end,:);
                u3=NaN(1440,720);
                u3(1:720,:)=u2;
                u3(721:end,:)=u1;
                % Narrow array down to the Gulf of Mexico
                lamin=nearestpoint(18, lat); lamax=nearestpoint(32,lat);
                lomin=nearestpoint(-98,lon); lomax=nearestpoint(-80,lon);
                plotlonSSS = lon(lomin:lomax); plotlatSSS = lat(lamin:lamax);
                [plotlonSSS,plotlatSSS] = meshgrid(plotlonSSS,plotlatSSS);
                plotSSS = u3(lomin:lomax,lamin:lamax);
                plotSSS = plotSSS';
                bigSSS(:,:,filenumSSS) = plotSSS;
                disp(filenumSSS);
                filenumSSS = filenumSSS + 1;
            end
        end
    else
        for j = 1:365
            lat = ncread(filenamesSSS2(j),'lat');
            lon = ncread(filenamesSSS2(j),'lon');
            sss = ncread(filenamesSSS2(j),'sss_smap');
            lon = double(lon); lat = double(lat); 
            % the lon is on 0 - 360 instead of -180 - 180
            dat1=lon(1:720);
            dat2=lon(721:end);
            dat3=dat2-360;
            dat4=NaN(1440,1);
            dat4(1:720)=dat3;
            dat4(721:end)=dat1;
            lon=dat4;
            % Rearrange Data to fit new longitudes
            u1=sss(1:720,:);
            u2=sss(721:end,:);
            u3=NaN(1440,720);
            u3(1:720,:)=u2;
            u3(721:end,:)=u1;
            % Narrow array down to the Gulf of Mexico
            lamin=nearestpoint(18, lat); lamax=nearestpoint(32,lat);
            lomin=nearestpoint(-98,lon); lomax=nearestpoint(-80,lon);
            plotlonSSS1 = lon(lomin:lomax); plotlatSSS1 = lat(lamin:lamax);
            [plotlonSSS,plotlatSSS] = meshgrid(plotlonSSS1,plotlatSSS1);
            plotSSS = u3(lomin:lomax,lamin:lamax);
            plotSSS = plotSSS';
            bigSSS(:,:,filenumSSS) = plotSSS;
            disp(filenumSSS);
            filenumSSS = filenumSSS + 1;
        end
    end
end
%% Saving SSS file
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/SSS_sat.mat';
save(filename_out,'bigSSS','plotlatSSS','plotlonSSS','plotlatSSS1','plotlonSSS1');
%% 
clear; clc;
cd /Volumes/Kate-Research/Data/Eddy_Extraction/Variables
load('SSS_sat.mat');
%% Compute anomalies
SSTsatref = NaN(57,73,365);
bigHYSST = bigSSS;
plotlatSST = plotlatSSS;
plotlonSST = plotlonSSS;
for i = 1:366
    meanholder = NaN(57,73,3);
    if i==60
        lpdayholder = NaN(57,73,2);
        lpdayholder(:,:,1) = bigHYSST(:,:,59);
        lpdayholder(:,:,2) = bigHYSST(:,:,60);
        mean19 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,1) = mean19;
        meanholder(:,:,2) = bigHYSST(:,:,425);
        lpdayholder = NaN(57,73,2);
        lpdayholder(:,:,1) = bigHYSST(:,:,790);
        lpdayholder(:,:,2) = bigHYSST(:,:,791);
        mean21 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,3) = mean21;
        lpdaymeanSST = mean(meanholder,3,'omitnan');
    elseif i<60
        meanholder(:,:,1) = bigHYSST(:,:,i);
        meanholder(:,:,2) = bigHYSST(:,:,(i+365));
        meanholder(:,:,3) = bigHYSST(:,:,(i+731));
        SSTsatref(:,:,i) = mean(meanholder,3,'omitnan');
    elseif i>60
        meanholder(:,:,1) = bigHYSST(:,:,i-1);
        meanholder(:,:,2) = bigHYSST(:,:,(i+365));
        meanholder(:,:,3) = bigHYSST(:,:,(i+730));
        SSTsatref(:,:,i-1) = mean(meanholder,3,'omitnan');
    end
end
SSTanom = NaN(57,73,1096);
for i = 1:1096
    if i == 425
        SSTanom(:,:,i) = bigHYSST(:,:,i)-lpdaymeanSST;
    elseif i<=365
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i);
    elseif i>=365 && i<425
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i-365);
    elseif i>425 && i<=731
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i-366);
    elseif i>=732 && i<=1096
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i-731);
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/SSSanom.mat';
save(filename_out,'SSTanom','-v7.3');
%% HYCOM Eddies
cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load('Eddy_Contours.mat');
% SSTanom = SSTanom;
AE_cont = AEs; CE_cont = CEs;
l = length(AE_cont(:,1,1));
SST_eddys_AE = cell(l,1096);
% Pull out SST for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    AE_cont_x = AE_cont(:,i,1); AE_cont_y = AE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        xa = double(AE_cont_x{j}); ya = double(AE_cont_y{j});
        SST = SSTanom(:,:,i);
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonSST,plotlatSST,xa,ya);% Take only datapoints within eddy boundary
%             outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
%             outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST and subtract the mean of the outside points to
            % get rid of seasonal variations
            SST_eddys_AE(j,i) = {SST};
        else 
            disp(j)
            disp(i)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(CE_cont(:,1,1));
SST_eddys_CE = cell(l,1096);
% Pull out SST for cyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    CE_cont_x = CE_cont(:,i,1); CE_cont_y = CE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        SST = SSTanom(:,:,i);
        xa = double(CE_cont_x{j}); ya = double(CE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonSST,plotlatSST,xa,ya);% Take only datapoints within eddy boundary
%             outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
%             outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST 
            SST_eddys_CE(j,i) = {SST};
        else 
            disp(j)
            disp(i)
        end
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SSSEddiesHYCOM.mat';
save(filename_out,'SST_eddys_AE','SST_eddys_CE','plotlonSST','plotlatSST','-v7.3');
clearvars -except SSTanom plotlonSST plotlatSST
% CMEMS Eddies
cd /Volumes/Kate-Research/MatlabPrograms/Eddy_Tracking_CMEMS/EXTRACTION/EDDY_PROPERTIES
load('Eddy_Contours.mat');
cd /Volumes/Kate-Research/Data/Eddy_Extraction
AE_cont = AEs; CE_cont = CEs;
l = length(AE_cont(:,1,1));
SST_eddys_AE = cell(l,1096);
% Pull out SST for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    AE_cont_x = AE_cont(:,i,1); AE_cont_y = AE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        SST = SSTanom(:,:,i);
        xa = double(AE_cont_x{j}); ya = double(AE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonSST,plotlatSST,xa,ya);% Take only datapoints within eddy boundary
            outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
            outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST 
            SST_eddys_AE(j,i) = {SST};
        else 
            disp(j)
            disp(i)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(CE_cont(:,1,1));
SST_eddys_CE = cell(l,1096);
% Pull out SST for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    CE_cont_x = CE_cont(:,i,1); CE_cont_y = CE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        SST = SSTanom(:,:,i);
        xa = double(CE_cont_x{j}); ya = double(CE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonSST,plotlatSST,xa,ya);% Take only datapoints within eddy boundary
            outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
            outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST 
            SST_eddys_CE(j,i) = {SST};
        else 
            disp(j)
            disp(i)
        end
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SSSEddiesCMEMS.mat';
save(filename_out,'SST_eddys_AE','SST_eddys_CE','plotlonSST','plotlatSST','-v7.3');
clear; clc;
%% Take time average of CMEMS/HYCOM for CEs and AEs
clear; clc;
cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load('SSSEddiesHYCOM.mat');
% HYCOM AEs
mean1 = NaN(21,1);
SSTbigmeanHYCOM_AE = NaN(21,1096);
for i = 1:1096
    AE = SST_eddys_AE(:,i);
    for j = 1:13
        if isempty(SST_eddys_AE(j,i))==0
            % take average SST from said eddy
            AE1 = double(AE{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1(j) = AE_mean;
        else
        end
    end
%     meannew = median(mean1,'omitnan');
    SSTbigmeanHYCOM_AE(:,i) = mean1;
end
% HYCOM CEs
mean1 = NaN(21,1);
SSTbigmeanHYCOM_CE = NaN(21,1096);
for i = 1:1096
    CE = SST_eddys_CE(:,i);
    for j = 1:18
        if isempty(SST_eddys_CE(j,i))==0
            % take average SST from said eddy
            CE1 = double(CE{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1(j) = CE_mean;
        else
        end
    end
%     meannew = median(mean1,'omitnan');
    SSTbigmeanHYCOM_CE(:,i) = mean1;
end
% Do the same for CMEMS
clearvars -except SSTbigmeanHYCOM_CE SSTbigmeanHYCOM_AE
renvar SSTbigmeanHYCOM_AE SSSbigmeanHYCOM_AE
renvar SSTbigmeanHYCOM_CE SSSbigmeanHYCOM_CE
load('SSSEddiesCMEMS.mat');
% HYCOM AEs
mean1 = NaN(21,1);
SSTbigmeanCMEMS_AE = NaN(21,1096);
for i = 1:1096
    AE = SST_eddys_AE(:,i);
    for j = 1:13
        if isempty(SST_eddys_AE(j,i))==0
            % take average SST from said eddy
            AE1 = double(AE{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1(j) = AE_mean;
        else
        end
    end
%     meannew = median(mean1,'omitnan');
    SSTbigmeanCMEMS_AE(:,i) = mean1;
end
% HYCOM CEs
mean1 = NaN(21,1);
SSTbigmeanCMEMS_CE = NaN(21,1096);
for i = 1:1096
    CE = SST_eddys_CE(:,i);
    for j = 1:18
        if isempty(SST_eddys_CE(j,i))==0
            % take average SST from said eddy
            CE1 = double(CE{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1(j) = CE_mean;
        else
        end
    end
%     meannew = median(mean1,'omitnan');
    SSTbigmeanCMEMS_CE(:,i) = mean1;
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SSS_eddies_avg.mat';
renvar SSTbigmeanCMEMS_CE SSSbigmeanCMEMS_CE
renvar SSTbigmeanCMEMS_AE SSSbigmeanCMEMS_AE
save(filename_out,'SSSbigmeanCMEMS_CE','SSSbigmeanCMEMS_AE','SSSbigmeanHYCOM_AE','SSSbigmeanHYCOM_CE','-v7.3');
%% PART TWO Model Data
% MODEL DATA!!!!!!!!!!!!! (SSS & SST)
clear variables; clc;
years = ["2019",'2020',"2021"];
pathtoHYCOMdata = '/Volumes/LACIE-GOM/GOM-HYCOM-1993-2020-Daliy3D-Ebenezer/daily'; % 1993-2020
pathtoHYCOMdata2021 = '/Volumes/Kate-Research/Data/HYCOM_2021/2021';
totalfilesHYCOM = 0;
for i = 1:2 % cycles through the years!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HYCOMpath = (pathtoHYCOMdata);
    yrstr=num2str(years(i));
    pathtoyearly = [HYCOMpath '/' yrstr '/'];
    cd(pathtoyearly);
    list = dir('*.nc');
    filenamesHYCOM = string({list.name});
    filenamesHYCOM = filenamesHYCOM';
    daysinfileHYCOM = length(filenamesHYCOM);
    totalfilesHYCOM = totalfilesHYCOM + daysinfileHYCOM;
end
cd(pathtoHYCOMdata2021);
list = dir('*.nc');
filenamesHYCOM = string({list.name});
filenamesHYCOM = filenamesHYCOM';
daysinfileHYCOM = length(filenamesHYCOM);
totalfilesHYCOM = totalfilesHYCOM + daysinfileHYCOM;
for i = 1:length(years) % cycles through the years!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yearnum = str2double(years(i));
    if yearnum == 2021
        cd(pathtoHYCOMdata2021);
        list = dir('*.nc');
        filenamesHYCOM = string({list.name});
        filenamesHYCOM = filenamesHYCOM';
        daysinfileHYCOM = length(filenamesHYCOM);
    else
        % Get HYCOM filenames for the specific year
        HYCOMpath = (pathtoHYCOMdata);
        yrstr=num2str(years(i));
        pathtoyearly = [HYCOMpath '/' yrstr '/'];
        cd(pathtoyearly);
        list = dir('*.nc');
        filenamesHYCOM = string({list.name});
        filenamesHYCOM = filenamesHYCOM';
        daysinfileHYCOM = length(filenamesHYCOM);
    end
    % Load in HYCOM data
    filenumHYCOM = 1;
    bigHYSST = NaN(346,451,totalfilesHYCOM);
    bigHYSSS = NaN(346,451,totalfilesHYCOM);
    plotHYSST= NaN(346,451);
    plotHYSSS = NaN(346,451);
    for j=1:daysinfileHYCOM
        % Load in HYCOM data
        if yearnum == 2021
            cd(pathtoHYCOMdata2021);
        else
            cd([HYCOMpath '/' yrstr '/']);
        end
        lat = ncread(filenamesHYCOM(j),'latitude'); % degrees north 90 to -90
        lon = ncread(filenamesHYCOM(j),'longitude'); % degrees east -180 to 180
        lon = double(lon); lat = double(lat);
        lamin=nearestpoint(18, lat); lamax=nearestpoint(32,lat);
        lomin=nearestpoint(-98,lon); lomax=nearestpoint(-80,lon);
        plotlonHYCOM = lon(lomin:lomax); plotlatHYCOM = lat(lamin:lamax);
        [plotlonHYCOM,plotlatHYCOM] = meshgrid(plotlonHYCOM,plotlatHYCOM);
        SST = ncread(filenamesHYCOM(j),'water_temp');
        SSS = ncread(filenamesHYCOM(j),'salinity');
        depth = ncread(filenamesHYCOM(j),'depth');
        depth = depth';
        plotHYSST = SST(lomin:lomax,lamin:lamax,1);
        plotHYSST = squeeze(plotHYSST);
        plotHYSST = plotHYSST';
        plotHYSSS = SSS(lomin:lomax,lamin:lamax,1);
        plotHYSSS = squeeze(plotHYSSS);
        plotHYSSS = plotHYSSS';
        % Create array to hold data to average later
        bigHYSST(:,:,filenumHYCOM) = plotHYSST;
        bigHYSSS(:,:,filenumHYCOM) = plotHYSSS;
        filenumHYCOM = filenumHYCOM + 1;
    end
end
%% Saving SST file
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/SST_SSS_HYCOM.mat';
save(filename_out,'bigHYSSS','bigHYSST','plotlatHYCOM','plotlonHYCOM');
clear; clc;
cd /Volumes/Kate-Research/Data/Eddy_Extraction/Variables
load('SST_SSS_HYCOM.mat');
%% Compute anomalies
SSTsatref = NaN(346,451,365);
SSSsatref = NaN(346,451,365);
for i = 1:366
    if i==60
        % SST
        meanholder = NaN(346,451,3);
        lpdayholder = NaN(346,451,2);
        lpdayholder(:,:,1) = bigHYSST(:,:,59);
        lpdayholder(:,:,2) = bigHYSST(:,:,60);
        mean19 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,1) = mean19;
        meanholder(:,:,2) = bigHYSST(:,:,425);
        lpdayholder = NaN(346,451,2);
        lpdayholder(:,:,1) = bigHYSST(:,:,790);
        lpdayholder(:,:,2) = bigHYSST(:,:,791);
        mean21 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,3) = mean21;
        lpdaymeanSST = mean(meanholder,3,'omitnan');
        % SSS
        meanholder = NaN(346,451,3);
        lpdayholder = NaN(346,451,2);
        lpdayholder(:,:,1) = bigHYSSS(:,:,59);
        lpdayholder(:,:,2) = bigHYSSS(:,:,60);
        mean19 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,1) = mean19;
        meanholder(:,:,2) = bigHYSSS(:,:,425);
        lpdayholder = NaN(346,451,2);
        lpdayholder(:,:,1) = bigHYSSS(:,:,790);
        lpdayholder(:,:,2) = bigHYSSS(:,:,791);
        mean21 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,3) = mean21;
        lpdaymeanSSS = mean(meanholder,3,'omitnan');
    elseif i<60
        % SST
        meanholder = NaN(346,451,3);
        meanholder(:,:,1) = bigHYSST(:,:,i);
        meanholder(:,:,2) = bigHYSST(:,:,(i+365));
        meanholder(:,:,3) = bigHYSST(:,:,(i+731));
        SSTsatref(:,:,i) = mean(meanholder,3,'omitnan');
        % SSS
        meanholder = NaN(346,451,3);
        meanholder(:,:,1) = bigHYSSS(:,:,i);
        meanholder(:,:,2) = bigHYSSS(:,:,(i+365));
        meanholder(:,:,3) = bigHYSSS(:,:,(i+731));
        SSSsatref(:,:,i) = mean(meanholder,3,'omitnan');
    elseif i>60
        % SST
        meanholder = NaN(346,451,3);
        meanholder(:,:,1) = bigHYSST(:,:,i-1);
        meanholder(:,:,2) = bigHYSST(:,:,(i+365));
        meanholder(:,:,3) = bigHYSST(:,:,(i+730));
        SSTsatref(:,:,i-1) = mean(meanholder,3,'omitnan');
        % SSS
        meanholder = NaN(346,451,3);
        meanholder(:,:,1) = bigHYSSS(:,:,i-1);
        meanholder(:,:,2) = bigHYSSS(:,:,(i+365));
        meanholder(:,:,3) = bigHYSSS(:,:,(i+730));
        SSSsatref(:,:,i-1) = mean(meanholder,3,'omitnan');
    end
end
%%
SSSanom = NaN(346,451,1096);
SSTanom = NaN(346,451,1096);
for i = 1:1096
    if i == 425
        SSTanom(:,:,i) = bigHYSST(:,:,i)-lpdaymeanSST;
        SSSanom(:,:,i) = bigHYSSS(:,:,i)-lpdaymeanSSS;
    elseif i<=365
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i);
        SSSanom(:,:,i) = bigHYSSS(:,:,i)-SSSsatref(:,:,i);
    elseif i>=365 && i<425
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i-365);
        SSSanom(:,:,i) = bigHYSSS(:,:,i)-SSSsatref(:,:,i-365);
    elseif i>425 && i<=731
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i-366);
        SSSanom(:,:,i) = bigHYSSS(:,:,i)-SSSsatref(:,:,i-366);
    elseif i>=732 && i<=1096
        SSTanom(:,:,i) = bigHYSST(:,:,i)-SSTsatref(:,:,i-731);
        SSSanom(:,:,i) = bigHYSSS(:,:,i)-SSSsatref(:,:,i-731);
    end
end
%% HYCOM Eddies
cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load('Eddy_Contours.mat');
AE_cont = AEs; CE_cont = CEs;
l = length(AE_cont(:,1,1));
SST_eddys_AE = cell(l,1096);
SSS_eddys_AE = cell(l,1096);
% Pull out SST for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    AE_cont_x = AE_cont(:,i,1); AE_cont_y = AE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        xa = double(AE_cont_x{j}); ya = double(AE_cont_y{j});
        SST = SSTanom(:,:,i);
        SSS = SSSanom(:,:,i);
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlatHYCOM,plotlonHYCOM,xa,ya);% Take only datapoints within eddy boundary
%             outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
            SSS(in==0) = NaN;
%             outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST and subtract the mean of the outside points to
            % get rid of seasonal variations
            SST_eddys_AE(j,i) = {SST};
            SSS_eddys_AE(j,i) = {SSS};
        else 
            disp(j)
            disp(i)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(CE_cont(:,1,1));
SST_eddys_CE = cell(l,1096);
SSS_eddys_CE = cell(l,1096);
% Pull out SST for cyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    CE_cont_x = CE_cont(:,i,1); CE_cont_y = CE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        SST = SSTanom(:,:,i);
        SSS = SSSanom(:,:,i);
        xa = double(CE_cont_x{j}); ya = double(CE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlatHYCOM,plotlonHYCOM,xa,ya);% Take only datapoints within eddy boundary
%             outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
            SSS(in==0) = NaN;
%             outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST 
            SST_eddys_CE(j,i) = {SST};
            SSS_eddys_CE(j,i) = {SSS};
        else 
            disp(j)
            disp(i)
        end
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SST_SSS_EddiesHYCOM.mat';
save(filename_out,'SST_eddys_AE','SSS_eddys_AE','SSS_eddys_CE','SST_eddys_CE','plotlatHYCOM','plotlonHYCOM','-v7.3');
% CMEMS Eddies
clearvars SST_eddys_AE SSS_eddys_AE SSS_eddys_CE SST_eddys_CE
cd /Volumes/Kate-Research/MatlabPrograms/Eddy_Tracking_CMEMS/EXTRACTION/EDDY_PROPERTIES
load('Eddy_Contours.mat');
cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load('SSTEddiesHYCOM.mat');
AE_cont = AEs; CE_cont = CEs;
l = length(AE_cont(:,1,1));
SST_eddys_AE = cell(l,1096);
SSS_eddys_AE = cell(l,1096);
% Pull out SST for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    AE_cont_x = AE_cont(:,i,1); AE_cont_y = AE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        SST = SSTanom(:,:,i);
        SSS = SSSanom(:,:,i);
        xa = double(AE_cont_x{j}); ya = double(AE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlatHYCOM,plotlonHYCOM,xa,ya);% Take only datapoints within eddy boundary
%             outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
            SSS(in==0) = NaN;
%             outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST 
            SST_eddys_AE(j,i) = {SST};
            SSS_eddys_AE(j,i) = {SSS};
        else 
            disp(j)
            disp(i)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(CE_cont(:,1,1));
SST_eddys_CE = cell(l,1096);
SSS_eddys_CE = cell(l,1096);
% Pull out SST for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    CE_cont_x = CE_cont(:,i,1); CE_cont_y = CE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        SST = SSTanom(:,:,i);
        SSS = SSSanom(:,:,i);
        xa = double(CE_cont_x{j}); ya = double(CE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonSST,plotlatSST,xa,ya);% Take only datapoints within eddy boundary
%             outside = SST; % create array to hold the SST data outside the eddy contour
            SST(in==0) = NaN; % sets all SST points outside eddy countour as NaNs!
            SSS(in==0) = NaN;
%             outside(in==1) = NaN; % sets all SST points inside the eddys to NaNs!
            % Save SST 
            SST_eddys_CE(j,i) = {SST};
            SSS_eddys_CE(j,i) = {SSS};
        else 
            disp(j)
            disp(i)
        end
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SST_SSS_EddiesCMEMS.mat';
save(filename_out,'SST_eddys_AE','SSS_eddys_AE','SSS_eddys_CE','SST_eddys_CE','plotlatHYCOM','plotlonHYCOM','-v7.3');
clear; clc;
%% Take time average of CMEMS/HYCOM for CEs and AEs
clear; clc;
cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load('SST_SSS_EddiesHYCOM.mat');
% HYCOM AEs
mean1sst = NaN(21,1);
mean1sss = NaN(21,1);
SSTbigmeanHYCOMx2_AE = NaN(21,1096);
SSSbigmeanHYCOMx2_AE = NaN(21,1096);
for i = 1:1096
    AEsst = SST_eddys_AE(:,i);
    AEsss = SSS_eddys_AE(:,i);
    for j = 1:13
        if isempty(SST_eddys_AE(j,i))==0
            % take average SST from said eddy
            AE1 = double(AEsst{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1sst(j) = AE_mean;
        else
        end
        if isempty(SSS_eddys_AE(j,i))==0
            % take average SST from said eddy
            AE1 = double(AEsss{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1sss(j) = AE_mean;
        else
        end
    end
    SSTbigmeanHYCOMx2_AE(:,i) = mean1sst;
    SSSbigmeanHYCOMx2_AE(:,i) = mean1sss;
end
% HYCOM CEs
mean1sst = NaN(21,1);
mean1sss = NaN(21,1);
SSTbigmeanHYCOMx2_CE = NaN(21,1096);
SSSbigmeanHYCOMx2_CE = NaN(21,1096);
for i = 1:1096
    CEsst = SST_eddys_CE(:,i);
    CEsss = SSS_eddys_CE(:,i);
    for j = 1:18
        if isempty(SST_eddys_CE(j,i))==0
            % take average SST from said eddy
            CE1 = double(CEsst{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1sst(j) = CE_mean;
        else
        end
        if isempty(SSS_eddys_CE(j,i))==0
            % take average SST from said eddy
            CE1 = double(CEsss{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1sss(j) = CE_mean;
        else
        end
    end
    SSTbigmeanHYCOMx2_CE(:,i) = mean1sst;
    SSSbigmeanHYCOMx2_CE(:,i) = mean1sss;
end
% Do the same for CMEMS
clearvars -except SSTbigmeanHYCOMx2_CE SSTbigmeanHYCOMx2_AE SSSbigmeanHYCOMx2_AE SSSbigmeanHYCOMx2_CE
load('SST_SSS_EddiesCMEMS.mat');
% HYCOM AEs
mean1sst = NaN(21,1);
mean1sss = NaN(21,1);
SSTHYbigmeanCMEMS_AE = NaN(21,1096);
SSSHYbigmeanCMEMS_AE = NaN(21,1096);
for i = 1:1096
    AEsst = SST_eddys_AE(:,i);
    AEsss = SSS_eddys_AE(:,i);
    for j = 1:13
        if isempty(SST_eddys_AE(j,i))==0
            % take average SST from said eddy
            AE1 = double(AEsst{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1sst(j) = AE_mean;
        else
        end
        if isempty(SSS_eddys_AE(j,i))==0
            % take average SST from said eddy
            AE1 = double(AEsss{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1sss(j) = AE_mean;
        else
        end
    end
    SSTHYbigmeanCMEMS_AE(:,i) = mean1sst;
    SSSHYbigmeanCMEMS_AE(:,i) = mean1sss;
end
% HYCOM CEs
mean1sst = NaN(21,1);
mean1sss = NaN(21,1);
SSTHYbigmeanCMEMS_CE = NaN(21,1096);
SSSHYbigmeanCMEMS_CE = NaN(21,1096);
for i = 1:1096
    CEsst = SST_eddys_CE(:,i);
    CEsss = SSS_eddys_CE(:,i);
    for j = 1:18
        if isempty(SST_eddys_CE(j,i))==0
            % take average SST from said eddy
            CE1 = double(CEsst{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1sst(j) = CE_mean;
        else
        end
        if isempty(SSS_eddys_CE(j,i))==0
            % take average SST from said eddy
            CE1 = double(CEsss{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1sss(j) = CE_mean;
        else
        end
    end
    SSTHYbigmeanCMEMS_CE(:,i) = mean1sst;
    SSSHYbigmeanCMEMS_CE(:,i) = mean1sss;
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/SST_SSS_eddies_avg_HYCOM.mat';
save(filename_out,'SSTHYbigmeanCMEMS_CE','SSTHYbigmeanCMEMS_AE','SSSHYbigmeanCMEMS_AE','SSSHYbigmeanCMEMS_CE','SSSbigmeanHYCOMx2_AE','SSSbigmeanHYCOMx2_CE','SSTbigmeanHYCOMx2_AE','SSTbigmeanHYCOMx2_CE','-v7.3');