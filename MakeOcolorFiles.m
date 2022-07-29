addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/Corinne_toolbox');
addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/m_map');
addpath('/Volumes/Lacie-SAN/SAN2/MATLAB/mytoolbox');% Sarahs toolbox :P
addpath('/Volumes/NEMO3/GOFFISH/Emily/ToolBox/cmocean_v2.0/cmocean');
addpath('/Volumes/Kate-Research/MatlabPrograms');
%% Make Ocolor files
clear; clc; %get rid of leftover spooky arrays :0
pathtoOcolor = '/Volumes/Lacie-SAN/SAN1/OceanColor/NOAA-DINEOF-gapfilled/'; % 2018-2022
pathtodata = '/Volumes/LACIE-GOM/GOM-HYCOM-1993-2020-Daliy3D-Ebenezer/daily/'; % 1993-2020
katepath = '/Volumes/Kate-Research/Data/Eddy_Extraction/';
output_dir= strcat(katepath);
years = ["2019",'2020',"2021"]; % List of years
%%
O_color_output = NaN(169,217,1096);
countfiles = 1;
for i=1:length(years)
    yearnum = str2double(years(i));
    filenamesO = [];
    yrstrO=num2str(years(i));
    filenamesO = [filenamesO;dir([pathtoOcolor yrstrO '/*.nc'])];
    % Create control date array
    datestart = datetime(yearnum,1,1);
    dateend = datetime(yearnum,12,31);
    controldate = datestart:dateend;
    controldate = controldate';
    controldays = length(controldate);
    skipcount = 0;
    for i=1:length(filenamesO)
        name =  filenamesO(i).name;
        odir = [pathtoOcolor yrstrO '/'];
        filenamesO1(i,:) = [odir name]; %% DT TO NRT!!!
    end
    renvar filenamesO1 filenamesO
    % Average missing days
    for i = 1:controldays
        if yearnum == 2020
            allmissing = [40;62;63;64;75;104;105;121;122;123;124;125;148;256;257];
            missing = [40;75;148];
            % need files 63, 104.5, 123, 121.5, 124.5, 256.5 (in order)
            filenames_O_61 = filenamesO(60,:);
            filenames_O_65 = filenamesO(61,:);
            Ocolor_before = double(ncread(filenames_O_61,'chlor_a'));
            Ocolor_after = double(ncread(filenames_O_65,'chlor_a'));
            O = NaN(4320,2160,2);
            O(:,:,1) = Ocolor_before;
            O(:,:,2) = Ocolor_after;
            Ocolor_63 = mean(O,3,'omitnan');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            filenames_O_103 = filenamesO(98,:);
            filenames_O_106 = filenamesO(99,:);
            Ocolor_before = double(ncread(filenames_O_103,'chlor_a'));
            Ocolor_after = double(ncread(filenames_O_106,'chlor_a')); 
            O = NaN(4320,2160,2);
            O(:,:,1) = Ocolor_before;
            O(:,:,2) = Ocolor_after;
            Ocolor_104_5 = mean(O,3,'omitnan');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            filenames_O_120 = filenamesO(113,:);
            filenames_O_126 = filenamesO(114,:);
            Ocolor_before = double(ncread(filenames_O_120,'chlor_a'));
            Ocolor_after = double(ncread(filenames_O_126,'chlor_a'));
            O = NaN(4320,2160,2);
            O(:,:,1) = Ocolor_before;
            O(:,:,2) = Ocolor_after;
            Ocolor_123 = mean(O,3,'omitnan');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Ocolor_before = double(ncread(filenames_O_120,'chlor_a'));
            Ocolor_after = Ocolor_123;
            O = NaN(4320,2160,2);
            O(:,:,1) = Ocolor_before;
            O(:,:,2) = Ocolor_after;
            Ocolor_121_5 = mean(O,3,'omitnan');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Ocolor_after = double(ncread(filenames_O_126,'chlor_a'));
            Ocolor_before = Ocolor_123;
            O = NaN(4320,2160,2);
            O(:,:,1) = Ocolor_before;
            O(:,:,2) = Ocolor_after;
            Ocolor_124_5 = mean(O,3,'omitnan');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            filenames_O_255 = filenamesO(242,:);
            filenames_O_258 = filenamesO(243,:);
            Ocolor_before = double(ncread(filenames_O_255,'chlor_a'));
            Ocolor_after = double(ncread(filenames_O_258,'chlor_a'));
            O = NaN(4320,2160,2);
            O(:,:,1) = Ocolor_before;
            O(:,:,2) = Ocolor_after;
            Ocolor_256_5 = mean(O,3,'omitnan');
            if ismember(i,missing) == 1
                filenames_O_before = filenamesO((i-skipcount-1),:);
                filenames_O_after = filenamesO(i-skipcount,:);
                Ocolor_before = double(ncread(filenames_O_before,'chlor_a'));
                Ocolor_after = double(ncread(filenames_O_after,'chlor_a'));
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==62
                Ocolor_before = double(ncread(filenames_O_61,'chlor_a'));
                Ocolor_after = Ocolor_63;
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==63
                Ocolor = Ocolor_63;
                skipcount = skipcount + 1;
            elseif i==64
                Ocolor_before = Ocolor_63;
                Ocolor_after = double(ncread(filenames_O_65,'chlor_a'));
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==104
                Ocolor_before = double(ncread(filenames_O_103,'chlor_a'));
                Ocolor_after = Ocolor_104_5;
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==105
                Ocolor_before = Ocolor_104_5;
                Ocolor_after = double(ncread(filenames_O_106,'chlor_a'));
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==121
                Ocolor_before = double(ncread(filenames_O_120,'chlor_a'));
                Ocolor_after = Ocolor_121_5;
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==122
                Ocolor_before = Ocolor_121_5;
                Ocolor_after = Ocolor_123;
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==123
                Ocolor = Ocolor_123;
                skipcount = skipcount + 1;
            elseif i==124
                Ocolor_before = Ocolor_123;
                Ocolor_after = Ocolor_124_5;
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==125
                Ocolor_before = Ocolor_124_5;
                Ocolor_after = double(ncread(filenames_O_126,'chlor_a'));
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==256
                Ocolor_before = double(ncread(filenames_O_255,'chlor_a'));
                Ocolor_after = Ocolor_256_5;
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif i==257
                Ocolor_before = Ocolor_256_5;
                Ocolor_after = double(ncread(filenames_O_258,'chlor_a'));
                O = NaN(4320,2160,2);
                O(:,:,1) = Ocolor_before;
                O(:,:,2) = Ocolor_after;
                Ocolor = mean(O,3,'omitnan');
                skipcount = skipcount + 1;
            elseif ismember(i,allmissing) == 0
                % Display file treated name and create associated filename
                filenames_O = filenamesO((i-skipcount),:);
                Ocolor = double(ncread(filenames_O,'chlor_a'));
            end
        else
            filenames_O = filenamesO(i,:);
            Ocolor = double(ncread(filenames_O,'chlor_a'));
        end
        Ocolorlat = double(ncread(filenames_O,'lat')); % degrees north 90 to -90
        Ocolorlon = double(ncread(filenames_O,'lon'));
        lamin=nearestpoint(18, Ocolorlat); lamax=nearestpoint(32,Ocolorlat);
        lomin=nearestpoint(-98,Ocolorlon); lomax=nearestpoint(-80,Ocolorlon);
        plotlonO = Ocolorlon(lomin:lomax); plotlatO = Ocolorlat(lamax:lamin);
        [plotlonO,plotlatO] = meshgrid(plotlonO,plotlatO);
        plotOcolor = Ocolor(lomin:lomax,lamax:lamin);
        plotOcolor = plotOcolor';
        O_color_output(:,:,countfiles) = plotOcolor;
        countfiles = countfiles+1;
        disp(filenames_O);
    end
end
% %% Masking 
% O_color_output1 = NaN(169,217,1096);
% test=O_color_output(:,:,560);
% figure
% pcolor(plotlonO,plotlatO,test);shading interp;hold on
% xlim([-100 -80]);ylim([16 35]);colorbar;colormap(jet(20));%caxis([15 38]);
% % take out weird stuff
% [x, y] = ginput;
% %[LONG_gin,LAT_gin]=m_xy2ll(x(:),y(:));
% %next
% y(end+1)=y(1);x(end+1)=x(1);
% [in,on]=inpolygon(plotlonO,plotlatO,x,y);% Take only datapoints within ORAS5 boundary
% for i = 1:1096
%     HYSSS = O_color_output(:,:,i);
%     HYSSS(in==0) = NaN;
%     O_color_output1(:,:,i) = HYSSS;
% end
%% Saving Ocolor files
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/Ocolor.mat';
save(filename_out,'O_color_output','plotlatO','plotlonO');
clear; clc;
%% Yay! I have Ocean color for 19-21. Now it's time to take the location of eddys from HYCOM/CMEMS and determine the ocean characteristics within them.
clear; clc;
cd /Volumes/Kate-Research/Data/Eddy_Extraction
load('Ocolor.mat');
ocolor = O_color_output;
clearvars O_color_output
%% Compute anomalies
ocolorref = NaN(169,217,365);
for i = 1:366
    meanholder = NaN(169,217,3);
    if i==60
        lpdayholder = NaN(169,217,2);
        lpdayholder(:,:,1) = ocolor(:,:,59);
        lpdayholder(:,:,2) = ocolor(:,:,60);
        mean19 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,1) = mean19;
        meanholder(:,:,2) = ocolor(:,:,425);
        lpdayholder = NaN(169,217,2);
        lpdayholder(:,:,1) = ocolor(:,:,790);
        lpdayholder(:,:,2) = ocolor(:,:,791);
        mean21 = mean(lpdayholder,3,'omitnan');
        meanholder(:,:,3) = mean21;
        lpdaymean = mean(meanholder,3,'omitnan');
    elseif i<60
        meanholder(:,:,1) = ocolor(:,:,i);
        meanholder(:,:,2) = ocolor(:,:,(i+365));
        meanholder(:,:,3) = ocolor(:,:,(i+731));
        ocolorref(:,:,i) = mean(meanholder,3,'omitnan');
    elseif i>60
        meanholder(:,:,1) = ocolor(:,:,i-1);
        meanholder(:,:,2) = ocolor(:,:,(i+365));
        meanholder(:,:,3) = ocolor(:,:,(i+730));
        ocolorref(:,:,i-1) = mean(meanholder,3,'omitnan');
    end
end
% Create control date array to add a lil structure
years = ["2019","2020","2021"];
datestart = datetime(2019,1,1);
dateend = datetime(2021,12,31);
controldate = datestart:dateend;
controldate = controldate';
controldays = length(controldate);
m = month(controldate); m = num2str(m);
d = day(controldate); d = num2str(d);
ocoloranom = NaN(169,217,1096);
for i = 1:1096
    if i == 425
        ocoloranom(:,:,i) = ocolor(:,:,i)-lpdaymean;
    elseif i<=365
        ocoloranom(:,:,i) = ocolor(:,:,i)-ocolorref(:,:,i);
    elseif i>=365 && i<425
        ocoloranom(:,:,i) = ocolor(:,:,i)-ocolorref(:,:,i-365);
    elseif i>425 && i<=731
        ocoloranom(:,:,i) = ocolor(:,:,i)-ocolorref(:,:,i-366);
    elseif i>=732 && i<=1096
        ocoloranom(:,:,i) = ocolor(:,:,i)-ocolorref(:,:,i-731);
    end
end
%% HYCOM Eddies
cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load('Eddy_Contours.mat');
O_color_output = ocoloranom;
AE_cont = AEs; CE_cont = CEs;
l = length(AE_cont(:,1,1));
Ocolor_eddys_AE = cell(l,1096);
% Pull out Ocolor for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    AE_cont_x = AE_cont(:,i,1); AE_cont_y = AE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        xa = double(AE_cont_x{j}); ya = double(AE_cont_y{j});
        Ocolor = O_color_output(:,:,i);
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonO,plotlatO,xa,ya);% Take only datapoints within eddy boundary
%             outside = Ocolor; % create array to hold the ocolor data outside the eddy contour
            Ocolor(in==0) = NaN; % sets all ocolor points outside eddy countour as NaNs!
%             outside(in==1) = NaN; % sets all ocolor points inside the eddys to NaNs!
            % Save Ocolor and subtract the mean of the outside points to
            % get rid of seasonal variations
            Ocolor_eddys_AE(j,i) = {Ocolor};
        else 
            disp(j)
            disp(i)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(CE_cont(:,1,1));
Ocolor_eddys_CE = cell(l,1096);
% Pull out Ocolor for cyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    CE_cont_x = CE_cont(:,i,1); CE_cont_y = CE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        Ocolor = O_color_output(:,:,i);
        xa = double(CE_cont_x{j}); ya = double(CE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonO,plotlatO,xa,ya);% Take only datapoints within eddy boundary
            outside = Ocolor; % create array to hold the ocolor data outside the eddy contour
            Ocolor(in==0) = NaN; % sets all ocolor points outside eddy countour as NaNs!
            outside(in==1) = NaN; % sets all ocolor points inside the eddys to NaNs!
            % Save Ocolor 
            Ocolor_eddys_CE(j,i) = {Ocolor};
        else 
            disp(j)
            disp(i)
        end
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/OcolorEddiesHYCOM.mat';
save(filename_out,'Ocolor_eddys_AE','Ocolor_eddys_CE','plotlonO','plotlatO','-v7.3');
% CMEMS Eddies
cd /Volumes/Kate-Research/MatlabPrograms/Eddy_Tracking_CMEMS/EXTRACTION/EDDY_PROPERTIES
load('Eddy_Contours.mat');
cd /Volumes/Kate-Research/Data/Eddy_Extraction
AE_cont = AEs(:,(8036:9131),:); CE_cont = CEs(:,(8036:9131),:);
l = length(AE_cont(:,1,1));
Ocolor_eddys_AE = cell(l,1096);
% Pull out Ocolor for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    AE_cont_x = AE_cont(:,i,1); AE_cont_y = AE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        Ocolor = O_color_output(:,:,i);
        xa = double(AE_cont_x{j}); ya = double(AE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonO,plotlatO,xa,ya);% Take only datapoints within eddy boundary
            outside = Ocolor; % create array to hold the ocolor data outside the eddy contour
            Ocolor(in==0) = NaN; % sets all ocolor points outside eddy countour as NaNs!
            outside(in==1) = NaN; % sets all ocolor points inside the eddys to NaNs!
            % Save Ocolor 
            Ocolor_eddys_AE(j,i) = {Ocolor};
        else 
            disp(j)
            disp(i)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l = length(CE_cont(:,1,1));
Ocolor_eddys_CE = cell(l,1096);
% Pull out Ocolor for Anticyclonic eddies
for i = 1:1096 % loop through days in 19-21
    % grabbing all the X and Y locations of AE and CE eddy contours
    % on that day
    CE_cont_x = CE_cont(:,i,1); CE_cont_y = CE_cont(:,i,2); 
    for j= 1:l % Loops through number of eddies on particular day
        Ocolor = O_color_output(:,:,i);
        xa = double(CE_cont_x{j}); ya = double(CE_cont_y{j});
        if isempty(xa)==0 && isempty(ya)==0
            % Pull only data points within the contour boundary
            [in,on]=inpolygon(plotlonO,plotlatO,xa,ya);% Take only datapoints within eddy boundary
            outside = Ocolor; % create array to hold the ocolor data outside the eddy contour
            Ocolor(in==0) = NaN; % sets all ocolor points outside eddy countour as NaNs!
            outside(in==1) = NaN; % sets all ocolor points inside the eddys to NaNs!
            % Save Ocolor 
            Ocolor_eddys_CE(j,i) = {Ocolor};
        else 
            disp(j)
            disp(i)
        end
    end
end
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/OcolorEddiesCMEMS.mat';
save(filename_out,'Ocolor_eddys_AE','Ocolor_eddys_CE','plotlonO','plotlatO','-v7.3');
clear; clc;
%% Take time average of CMEMS/HYCOM for CEs and AEs
clear; clc;
cd /Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES
load('OcolorEddiesHYCOM.mat');
% HYCOM AEs
mean1 = NaN(13,1);
bigmeanHYCOM_AE = NaN(1096,1);
for i = 1:1096
    AE = Ocolor_eddys_AE(:,i);
    for j = 1:13
        if isempty(Ocolor_eddys_AE(j,i))==0
            % take average Ocolor from said eddy
            AE1 = double(AE{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1(j) = AE_mean;
        else
        end
    end
    meannew = median(mean1,'omitnan');
    bigmeanHYCOM_AE(i) = meannew;
end
% HYCOM CEs
mean1 = NaN(18,1);
bigmeanHYCOM_CE = NaN(1096,1);
for i = 1:1096
    CE = Ocolor_eddys_CE(:,i);
    for j = 1:18
        if isempty(Ocolor_eddys_CE(j,i))==0
            % take average Ocolor from said eddy
            CE1 = double(CE{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1(j) = CE_mean;
        else
        end
    end
    meannew = median(mean1,'omitnan');
    bigmeanHYCOM_CE(i) = meannew;
end
% Do the same for CMEMS
clearvars -except bigmeanHYCOM_CE bigmeanHYCOM_AE
load('OcolorEddiesCMEMS.mat');
% HYCOM AEs
mean1 = NaN(13,1);
bigmeanCMEMS_AE = NaN(1096,1);
for i = 1:1096
    AE = Ocolor_eddys_AE(:,i);
    for j = 1:13
        if isempty(Ocolor_eddys_AE(j,i))==0
            % take average Ocolor from said eddy
            AE1 = double(AE{j});
            AE_mean = mean(AE1,'all','omitnan');
            mean1(j) = AE_mean;
        else
        end
    end
    meannew = median(mean1,'omitnan');
    bigmeanCMEMS_AE(i) = meannew;
end
% HYCOM CEs
mean1 = NaN(18,1);
bigmeanCMEMS_CE = NaN(1096,1);
for i = 1:1096
    CE = Ocolor_eddys_CE(:,i);
    for j = 1:18
        if isempty(Ocolor_eddys_CE(j,i))==0
            % take average Ocolor from said eddy
            CE1 = double(CE{j});
            CE_mean = mean(CE1,'all','omitnan');
            mean1(j) = CE_mean;
        else
        end
    end
    meannew = median(mean1,'omitnan');
    bigmeanCMEMS_CE(i) = meannew;
end
% climanom_CMEMS_CE = mean(bigmeanCMEMS_CE,'all','omitnan');
% climanom_CMEMS_AE = mean(bigmeanCMEMS_AE,'all','omitnan');
% climanom_HYCOM_CE = mean(bigmeanHYCOM_CE,'all','omitnan');
% climanom_HYCOM_AE = mean(bigmeanHYCOM_AE,'all','omitnan');
filename_out = '/Volumes/Kate-Research/Data/Eddy_Extraction/EDDY_PROPERTIES/Ocolor_eddies_avg.mat';
save(filename_out,'bigmeanCMEMS_CE','bigmeanCMEMS_AE','bigmeanHYCOM_AE','bigmeanHYCOM_CE','-v7.3');