% For Mid Term Review
% Results section 3.1 and 3.2
% Analysis for the entire Southern Ocean
clear;
warning('off','all');

%% Process the data from PISCES
% Data files 
exp0ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\BAIT_test_1m_ptrc_Y2049.nc";
exp0ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\BAIT_test_1m_ptrc_Y2050.nc";
exp1ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\Si8Fe0.5_1m_ptrc_Y2599.nc";
exp1ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\Si8Fe0.5_1m_ptrc_Y2600.nc";
exp2ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\nFe0.3_1m_ptrc_Y2599.nc";
exp2ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\nFe0.3_1m_ptrc_Y2600.nc";
exp3ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\pFe0.1_1m_ptrc_Y2599.nc";
exp3ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\pFe0.1_1m_ptrc_Y2600.nc";
exp4ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\pFe0.1nFe0.3dFe0.5_1m_ptrc_Y2599.nc";
exp4ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\pFe0.1nFe0.3dFe0.5_1m_ptrc_Y2600.nc";
exp5ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\PIa10_1m_ptrc_Y2599.nc";
exp5ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\PIa10_1m_ptrc_Y2600.nc";
exp6ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\dia.aPI_1m_ptrc_Y2599.nc";
exp6ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\dia.aPI_1m_ptrc_Y2600.nc";

domain = "D:\OneDrive - Australian National University\PhD\Project1\NutrientLightColimitation\data\ORCA_R2_zps_domcfg_4.0.nc";

% Read data
lon_read = ncread(domain,"nav_lon");
lat_read = ncread(domain,"nav_lat");
depth_read = ncread(domain,"nav_lev");

lon = lon_read(:,1:50);
lat = lat_read(:,1:50);
depth = depth_read(1:19);

iceesp = 0.1;

mon = (0.5:11.5)';
monlist = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

explist = {'Control','Exp1: D-','Exp2: N-','Exp3: P-','Exp4: All-','Exp5: PI+','Exp6: D- PI+'};

ptrcvar = {'NCHL','DCHL','PCHL'};

% Data reading and pre-analysis
for i = 1:length(explist)
    ptrcfile1 = eval(['exp' num2str(i-1) 'ptrc1']);
    ptrcfile2 = eval(['exp' num2str(i-1) 'ptrc2']);
    
    for j = 1:length(ptrcvar)
        data_read1 = ncread(ptrcfile1,ptrcvar{j});
        data_read2 = ncread(ptrcfile2,ptrcvar{j});
        data1_nosi = seaicerm(data_read1,iceesp);
        data2_nosi = seaicerm(data_read2,iceesp);
        data_comb = shyear(data1_nosi(:,1:50,:,:),data2_nosi(:,1:50,:,:));
        eval([ptrcvar{j} '_exp' num2str(i-1) ' = data_comb;']);
    end
end

for i = 1:length(explist)
    eval(['TCHL_exp' num2str(i-1) '=NCHL_exp' num2str(i-1) '+DCHL_exp' num2str(i-1) '+PCHL_exp' num2str(i-1) ';']);
end

% DCM analysis
for i = 1:length(explist)
    data = eval(['TCHL_exp' num2str(i-1)]);
    [CHLM,CHLMind] = max(data,[],3); % Calculate the max
    CHLM(CHLM<=1.2*data(:,:,1,:)) = NaN;
    CHLM = squeeze(CHLM);
    CHLMind = squeeze(CHLMind);
    % chlorophyll concentration at each point and the depth index 
    DCM = CHLM; % Calculate the deep chlorophyll maxima, in other words, 
    DCMind = CHLMind;
    DCMind(isnan(DCM)) = NaN;
    DCMd = NaN(size(DCM));

    % Calculate DCMs depth 
    for ji = 1:size(DCMind,1)
        for jj = 1:size(DCMind,2)
            for jk = 1:size(DCMind,3)
                if isnan(DCMind(ji,jj,jk))
                    DCMd(ji,jj,jk) = NaN;
                else
                    DCMd(ji,jj,jk) = depth_read(DCMind(ji,jj,jk));
                end
            end
        end
    end

    % Calculate DCM month number
    DCMnumbers = squeeze(sum(~isnan(DCM),3));

    DCMmonper = DCMnumbers/12*100;

    % Count DCM pixels
    Count_DCM = squeeze(sum(sum(~isnan(DCM),2),1));

    tot_chl = squeeze(sum(sum(~isnan(data(:,:,1,:)),2),1));

    DCM_per_ssl = Count_DCM./tot_chl*100;

    % Calculate latitudinal DCM proportion variation
    DCM_per_lat = squeeze(sum(sum(~isnan(DCM),3),1))./squeeze(sum(sum(~isnan(data(:,:,1,:)),4),1))*100;

    % Export the variables
    eval(['DCM_exp' num2str(i-1) '=DCM;']);
    eval(['DCMind_exp' num2str(i-1) '=DCMind;']);
    eval(['DCMd_exp' num2str(i-1) '=DCMd;']);
    eval(['DCMmonper_exp' num2str(i-1) '=DCMmonper;']);
    eval(['DCMpts_exp' num2str(i-1) '=Count_DCM;']);
    eval(['DCM_per_ssl_exp' num2str(i-1) '=DCM_per_ssl;'])
    eval(['DCM_per_lat_exp' num2str(i-1) '=DCM_per_lat;'])
end

%% Process the data from Copernicus
% Read the grid
filename = ['D:\OneDrive - Australian National University\Data\copernicus\' ...
    '2020jul-2021jun\cmems_obs_glo_bgc3d_rep_weekly_20200701T0000Z_P20220621.nc'];
latcop_read = ncread(filename,'latitude');
loncop_read = ncread(filename,'longitude');
depcop_read = ncread(filename,'depth');

latcop = latcop_read(1:209);
loncop = loncop_read([1:1440 1]);

% Read chlorophyll data from files 
folderPath = 'D:\OneDrive - Australian National University\Data\copernicus\2020jul-2021jun';

filePattern = fullfile(folderPath,'*.nc');
ncFiles = dir(filePattern);

data = NaN(1441,209,36,52);
timeraw = NaN(length(ncFiles),1);

for k = 1:length(ncFiles)
    baseFileName = ncFiles(k).name;
    fullFileName = fullfile(ncFiles(k).folder,baseFileName);

    chlData = ncread(fullFileName,'chl');
    timeraw(k,1) = ncread(fullFileName,'time');

    data(:,:,:,k) = chlData([1:1440 1],1:209,:);
end

time = datetime(1950,1,1,0,0,0) + hours(timeraw);
dateNums = datenum(time);
dateVec = datevec(dateNums);
daysInMonth = eomday(dateVec(:,1), dateVec(:,2)); 
dayFraction = dateVec(:,3) ./ daysInMonth;
dateVec(dateVec(:,2)<7,2) = dateVec(dateVec(:,2)<7,2) + 12;
mon_cop = mon(dateVec(:,2) - 6) + dayFraction - 0.5; 

% Calculate the DCM occurance
data(data==0) = NaN;
[CHLM,CHLMind] = max(data,[],3); % Calculate the max
CHLM(CHLM<=1.2*data(:,:,1,:)) = NaN;
CHLM = squeeze(CHLM);
CHLMind = squeeze(CHLMind);
% chlorophyll concentration at each point and the depth index 
DCM = CHLM; % Calculate the deep chlorophyll maxima, in other words, 
DCMind = CHLMind;
DCMind(isnan(DCM)) = NaN;
DCMd = NaN(size(DCM));
% Calculate DCMs depth 
for ki = 1:size(DCMind,1)
    for kj = 1:size(DCMind,2)
        for kk = 1:size(DCMind,3)
            if isnan(DCMind(ki,kj,kk))
                DCMd(ki,kj,kk) = NaN;
            else
                DCMd(ki,kj,kk) = depcop_read(DCMind(ki,kj,kk));
            end
        end
    end
end
        
% Calculate DCM month number
DCMm = squeeze(sum(~isnan(DCM),3));
DCMm(DCMm == 0) = NaN;
DCMmonper = DCMm/53*100;

% Count DCM pixels
DCMpoints = squeeze(sum(sum(~isnan(DCM),2),1));
tot_chl = squeeze(sum(sum(~isnan(data(:,:,1,:)),2),1));
DCM_per_ssl = DCMpoints./tot_chl*100;
    
% Calculate latitudinal DCM proportion variation
DCM_per_lat = squeeze(sum(sum(~isnan(DCM),3),1))./squeeze(sum(sum(~isnan(data(:,:,1,:)),4),1))*100;

DCM_cop = DCM;
DCMind_cop = DCMind;
DCMd_cop = DCMd;
DCMmonper_cop = DCMmonper;
DCMpts_cop = DCMpoints;
DCM_per_ssl_cop = DCM_per_ssl;
DCM_per_lat_cop = DCM_per_lat;

[longrid_cop,latgrid_cop] = ndgrid(loncop,latcop);

% Make Plots
%% Results Fig 3.1.1 DCM map for all experiments including Copernicus
% For value check purpose
data1 = squeeze(DCM_exp0(:,:,7));
data2 = squeeze(DCM_exp1(:,:,7));
data3 = squeeze(DCM_exp2(:,:,7));
data4 = squeeze(DCM_exp3(:,:,7));
data5 = squeeze(DCM_exp4(:,:,7));
data6 = squeeze(DCM_exp5(:,:,7));
data7 = squeeze(DCM_exp6(:,:,7));
data8 = squeeze(DCM_cop(:,:,29));

min([data1(:);data2(:);data3(:);data4(:);data5(:);data6(:);data7(:);data8(:)]);
max([data1(:);data2(:);data3(:);data4(:);data5(:);data6(:);data7(:);data8(:)]);

plotnamelist = {'a) Control','b) Exp1: D Fe-','c) Exp2: N Fe-','d) Exp3: P Fe-',...
    'e) Exp4: DNP Fe-','f) Exp5: PI+','g) Exp6: D Fe- PI+','h) Copernicus'};
plotabbrlist = {'exp0','exp1','exp2','exp3','exp4','exp5','exp6','cop'};
plottimelist = [7;7;7;7;7;7;7;29];

for i = 1:length(plotnamelist)
    data = squeeze(eval(['DCM_' plotabbrlist{i} '(:,:,' num2str(plottimelist(i)) ')']));
    data(data>0.01 & data<=0.1) = (log10(data(data>0.01 & data<=0.1))+1)*0.2;
    data(data>0.1 & data<=0.2) = (data(data>0.1 & data<=0.2)-0.1)*2;
    data(data>1 & data<=1.5) = (data(data>1 & data<=1.5)-1)*0.4+1;
    data(data>1.5 & data<=2) = (data(data>1.5 & data<=2)-1.5)*0.4+1.2;
    data(data>2 & data<=5) = (data(data>2 & data<=5)-2)*0.0667+1.4;
    data(data>5 & data<=10) = (data(data>5 & data<=10)-1)*0.04+1.6;
    eval(['data' num2str(i) '=data;']);
end

splist = [1;2;4;5;7;8;10;11];

figure(1)
set(gcf,'Position',[100 50 1000 1600])
for i = 1:length(explist)
    data = eval(['data' num2str(i)]);
    subplot(4,3,splist(i))
    m_proj('stereographic','lat',-90,'long',0,'radius',60);
    m_pcolor(lon,lat,data)
    shading flat;
    colormap(chl_colmap);
    m_coast('patch',[.7 .7 .7]);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',6,'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle');
    title(plotnamelist{i})
    set(get(gca,'title'),'Position',[0 1.35 0])
    clim([-0.2 1.8])
    
end

% Subplot for Copernicus
i=8;
data = eval(['data' num2str(i)]);
subplot(4,3,splist(i))
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(longrid_cop,latgrid_cop,data)
shading flat;
colormap(chl_colmap);
m_coast('patch',[.7 .7 .7]);
m_coast('linewidth',1,'color','k');
m_grid('xtick',6,'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle');
title(plotnamelist{i})
set(get(gca,'title'),'Position',[0 1.35 0])
clim([-0.2 1.8])
c=colorbar;
c.Ticks = -0.2:0.2:1.8;
c.TickLabels = [0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2, 5, 10]; 
c.Label.String = 'Chlorophyll (\mug/L)';
set(c,'Location','eastoutside','Position',[0.7, 0.3, 0.02, 0.4])

exportgraphics(gcf,'D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Figure\fig3.1DCM_jan.jpg',...
    'Resolution',600);

%% Results Fig 3.1.2 DCM depth map for all experiments including Copernicus
for i = 1:length(plotnamelist)
    data = squeeze(eval(['DCMd_' plotabbrlist{i} '(:,:,' num2str(plottimelist(i)) ')']));
    eval(['data' num2str(i) '=data;']);
end

figure(2)
set(gcf,'Position',[100 50 1000 1600])
for i = 1:length(explist)
    data = eval(['data' num2str(i)]);
    subplot(4,3,splist(i))
    m_proj('stereographic','lat',-90,'long',0,'radius',60);
    m_pcolor(lon,lat,data)
    shading flat;
    colormap(flipud(m_colmap('blue',256)));
    m_coast('patch',[.7 .7 .7]);
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',6,'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle');
    title(plotnamelist{i})
    set(get(gca,'title'),'Position',[0 1.35 0])
    clim([0 200])
    
end

% Subplot for Copernicus
i=8;
data = eval(['data' num2str(i)]);
subplot(4,3,splist(i))
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(longrid_cop,latgrid_cop,data)
shading flat;
colormap(flipud(m_colmap('blue',256)));
m_coast('patch',[.7 .7 .7]);
m_coast('linewidth',1,'color','k');
m_grid('xtick',6,'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle');
title(plotnamelist{i})
set(get(gca,'title'),'Position',[0 1.35 0])
clim([0 200])
c=colorbar;
c.Label.String = 'Depth (m)';
set(c,'Location','eastoutside','Position',[0.7, 0.3, 0.02, 0.4])

exportgraphics(gcf,'D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Figure\fig3.2DCMd_jan.jpg',...
    'Resolution',600);

%% Results Fig 3.2.1 seasonal DCM proportion
plotnamelist = {'Control','Exp1: D Fe-','Exp2: N Fe-','Exp3: P Fe-',...
    'Exp4: DNP Fe-','Exp5: PI+','Exp6: D Fe- PI+','Copernicus'};

figure(3)
set(gcf,'Position',[100 100 1000 500])
plot(mon,DCM_per_ssl_exp0,'k','LineWidth',2)
hold on;
plot(mon,DCM_per_ssl_exp1,'g','LineWidth',2)
plot(mon,DCM_per_ssl_exp2,'b','LineWidth',2)
plot(mon,DCM_per_ssl_exp3,'c','LineWidth',2)
plot(mon,DCM_per_ssl_exp4,'r','LineWidth',2)
plot(mon,DCM_per_ssl_exp5,'y','LineWidth',2)
plot(mon,DCM_per_ssl_exp6,'m','LineWidth',2)
plot(mon_cop,DCM_per_ssl_cop,'k--','LineWidth',2)
xlim([0 12])
xticks(0:11)
xticklabels(monlist)
xlabel('Month')
ylabel('Grid points (%)')
set(gca,'FontSize',12)
legend(plotnamelist,'FontSize',12)
grid on;
exportgraphics(gcf,'D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Figure\fig3.3DCMprop_ssn.jpg',...
    'Resolution',600);

%% Results Fig 3.2.2 latitudinal DCM proportion
plotnamelist = {'Control','Exp1: D Fe-','Exp2: N Fe-','Exp3: P Fe-',...
    'Exp4: DNP Fe-','Exp5: PI+','Exp6: D Fe- PI+','Copernicus'};

figure(4)
set(gcf,'Position',[100 100 1000 500])
plot(abs(lat(1,:)),DCM_per_lat_exp0,'k','LineWidth',2)
hold on;
plot(abs(lat(1,:)),DCM_per_lat_exp1,'g','LineWidth',2)
plot(abs(lat(1,:)),DCM_per_lat_exp2,'b','LineWidth',2)
plot(abs(lat(1,:)),DCM_per_lat_exp3,'c','LineWidth',2)
plot(abs(lat(1,:)),DCM_per_lat_exp4,'r','LineWidth',2)
plot(abs(lat(1,:)),DCM_per_lat_exp5,'y','LineWidth',2)
plot(abs(lat(1,:)),DCM_per_lat_exp6,'m','LineWidth',2)
plot(abs(latgrid_cop(1,:)),DCM_per_lat_cop,'k--','LineWidth',2)
set(gca,"XDir","reverse")
xlabel('Latitude')
ylabel('Grid points (%)')
set(gca,'FontSize',12)
legend(plotnamelist,'FontSize',12,'Location','southeast')
grid on;
exportgraphics(gcf,'D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Figure\fig3.4DCMprop_lat.jpg',...
    'Resolution',600);


function output=shyear(lastyr,nextyr)
% combine two years datasets and keep July to December of last year and
% January to June of next year
output = lastyr;
output(:,:,:,1:6) = lastyr(:,:,:,7:12);
output(:,:,:,7:12) = nextyr(:,:,:,1:6);
end