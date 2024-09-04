% Mid-Term figures
% Southern Ocean fronts
clear;
warning('off','all');

% Data file
chlfile = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\BAIT_test_1m_ptrc_Y2050.nc";
domain = "D:\OneDrive - Australian National University\PhD\Project1\NutrientLightColimitation\data\ORCA_R2_zps_domcfg_4.0.nc";

eps=0.1;

lon_read = ncread(domain,"nav_lon");
lat_read = ncread(domain,"nav_lat");
nchl_read = ncread(chlfile,"NCHL");
dchl_read = ncread(chlfile,"DCHL");
pchl_read = ncread(chlfile,"PCHL");
chl_read = seaicerm(nchl_read + dchl_read + pchl_read,eps);

lon = lon_read(:,1:50);
lat = lat_read(:,1:50);
chl = squeeze(chl_read(:,1:50,1,1));

% Make customised chl for plotting
chl(chl>0.01 & chl<=0.1) = (log10(chl(chl>0.01 & chl<=0.1))+1)*0.2;
chl(chl>0.1 & chl<=0.2) = (chl(chl>0.1 & chl<=0.2)-0.1)*2;
chl(chl>1 & chl<=1.5) = (chl(chl>1 & chl<=1.5)-1)*0.4+1;
chl(chl>1.5 & chl<=2) = (chl(chl>1.5 & chl<=2)-1.5)*0.4+1.2;
chl(chl>2 & chl<=5) = (chl(chl>2 & chl<=5)-2)*0.0667+1.4;
chl(chl>5 & chl<=10) = (chl(chl>5 & chl<=10)-1)*0.04+1.6;

lon_saz = [135:140 140:-1:135 135];
lat_saz = [-47*ones(1,6) -51*ones(1,6) -47];

lon_pfz = [135:140 140:-1:135 135];
lat_pfz = [-53*ones(1,6) -55*ones(1,6) -53];

lon_az = [135:140 140:-1:135 135];
lat_az = [-58*ones(1,6) -61*ones(1,6) -58];

% Read front data
frontname = {'stf','saf','pf'};
for i = 1:length(frontname)
    txtfile = ['D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Data\' frontname{i} '.txt'];
    fid = fopen(txtfile,'r');
    header = fgetl(fid);
    data = textscan(fid,'%f %f');
    fclose(fid);
    line=[];
    line(:,1) = data{1};
    line(:,2) = data{2};
    eval(['line_' frontname{i} '=line;']);
end

figure(1);
set(gcf,'Position',[100 100 650 650])
subplot(1,1,1)
set(gca,'Position',[0.1 0.2 0.7 0.7])
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(lon,lat,chl)
shading flat;
colormap(chl_colmap);
m_grid('xtick',6,'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle');
m_line(line_stf(:,1),line_stf(:,2),'color','y','linewi',2)
m_line(line_saf(:,1),line_saf(:,2),'color','r','linewi',2)
m_line(line_pf(:,1),line_pf(:,2),'color','m','linewi',2)
m_line(lon_saz,lat_saz,'color','k','linewi',2)
m_line(lon_pfz,lat_pfz,'color','k','linewi',2)
m_line(lon_az,lat_az,'color','k','linewi',2)
m_coast('patch',[.7 .7 .7]);
m_coast('linewidth',1,'color','k');
clim([-0.2 1.8])
c=colorbar;
c.Ticks = -0.2:0.2:1.8;
c.TickLabels = [0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2, 5, 10]; 
c.Label.String = 'Chlorophyll (\mug/L)';
set(c,'Location','southoutside','Position',[0.15, 0.12, 0.6, 0.02])
text(-0.1,-0.82,'STF','FontWeight','bold')
text(-0.1,-0.59,'SAF','FontWeight','bold')
text(-0.1,-0.5,' PF','FontWeight','bold')

exportgraphics(gcf,'D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Figure\front.jpg','Resolution',600)