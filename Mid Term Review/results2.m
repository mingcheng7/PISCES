% For Mid Term Review
% Results section 3.3 and discussion
% Analysis for different zones
clear;
warning('off','all');

%% Process the data from PISCES
% Data files 
exp0ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\BAIT_test_1m_ptrc_Y2049.nc";
exp0ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\BAIT_test_1m_ptrc_Y2050.nc";
exp1ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\Si8Fe0.5_1m_ptrc_Y2599.nc";
exp1ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\Si8Fe0.5_1m_ptrc_Y2600.nc";
exp2ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\PIa10_1m_ptrc_Y2599.nc";
exp2ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\PIa10_1m_ptrc_Y2600.nc";
exp3ptrc1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\dia.aPI_1m_ptrc_Y2599.nc";
exp3ptrc2 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\dia.aPI_1m_ptrc_Y2600.nc";
exp0diad1 = "D:\OneDrive - Australian National University\PhD\Project1\SiFeSensitivity\Data\experiment_phyto\BAIT_test_1m_diad_Y2050.nc";

domain = "D:\OneDrive - Australian National University\PhD\Project1\NutrientLightColimitation\data\ORCA_R2_zps_domcfg_4.0.nc";

% Read data
lon_read = ncread(domain,"nav_lon");
lat_read = ncread(domain,"nav_lat");
depth_read = ncread(domain,"nav_lev");

lon = lon_read(:,1:50);
lat = lat_read(:,1:50);
depth = depth_read(1:19);

mon = (1:12)';
monlist = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

explist = {'Control','Exp1: D-','Exp5: PI+','Exp6: D- PI+'};

ptrcvar = {'NCHL','DCHL','PCHL','Si','Fer','NH4','PHY','PHY2','PIC'};

% Data reading and pre-analysis
for i = 1:length(explist)
    ptrcfile1 = eval(['exp' num2str(i-1) 'ptrc1']);
    ptrcfile2 = eval(['exp' num2str(i-1) 'ptrc2']);
    
    for j = 1:length(ptrcvar)
        data_read1 = ncread(ptrcfile1,ptrcvar{j});
        data_read2 = ncread(ptrcfile2,ptrcvar{j});
        data_comb = shyear(data_read1(:,1:50,:,:),data_read2(:,1:50,:,:));
        eval([ptrcvar{j} '_exp' num2str(i-1) ' = data_comb;']);
    end
end

MLD_read = ncread(exp0diad1,'MLD');
MLD = shyear2(MLD_read(:,1:50,:),MLD_read(:,1:50,:));

for i = 1:length(explist)
    eval(['TCHL_exp' num2str(i-1) '=NCHL_exp' num2str(i-1) '+DCHL_exp' num2str(i-1) '+PCHL_exp' num2str(i-1) ';']);
    eval(['TBioC_exp' num2str(i-1) '=PHY_exp' num2str(i-1) '+PHY2_exp' num2str(i-1) '+PIC_exp' num2str(i-1) ';']);
    eval(['rCChl_exp' num2str(i-1) '=TBioC_exp' num2str(i-1) './TCHL_exp' num2str(i-1) '*12;']);
end

% Extract these data in specific regions
allvar = {'TCHL','Fer','Si','NH4','rCChl'};

for i = 1:length(allvar)
    for j = 1:length(explist)
        data = eval([allvar{i} '_exp' num2str(j-1)]);
        data_saz = squeeze(mean(mean(data(30:32,37:39,1:19,:),2),1));
        data_pfz = squeeze(mean(mean(data(30:32,33:34,1:19,:),2),1));
        data_az = squeeze(mean(mean(data(30:32,28:30,1:19,:),2),1));
        eval([allvar{i} '_exp' num2str(j-1) '_saz=data_saz;']);
        eval([allvar{i} '_exp' num2str(j-1) '_pfz=data_pfz;']);
        eval([allvar{i} '_exp' num2str(j-1) '_az=data_az;']);
    end
end

MLD_saz = squeeze(mean(mean(MLD(30:32,37:39,:),2),1));
MLD_pfz = squeeze(mean(mean(MLD(30:32,33:34,:),2),1));
MLD_az = squeeze(mean(mean(MLD(30:32,28:30,:),2),1));

% Make plots for Section 3.3
%% Fig 3.3.1 Seasonal Chlorophyll distribution in Antarctic Zone
plotnamelist = {'a) Control','b) Exp1: D Fe-','c) Exp5: PI+','d) Exp6: D Fe- PI+'};

figure(1)
set(gcf,'Position',[200 50 2000 1000])
for i = 1:length(explist)
    data = eval(['TCHL_exp' num2str(i-1) '_az']);
    subplot(4,2,i*2-1)
    contourf(mon,depth,data,256,'linestyle','none')
    shading flat;
    colormap(flipud(m_colmap('green',256)));
    hold on;
    plot(mon,MLD_saz,'Color',[.3 .3 .3],'LineStyle','--','LineWidth',0.5)
    set(gca,'YDir','reverse')
    ylabel('Depth (m)')
    xlabel('Month')
    xticks(mon)
    xticklabels(monlist)
    clim([0 1.2])
    ylim([-inf 300])
    yticks(50:50:300)
    title(plotnamelist{i})
    if i == 1
        c=colorbar;
        c.Ticks = 0:0.2:1.2;
        c.Label.String = 'Chlorophyll (\mug/L)';
        set(c,'Location','eastoutside','Position',[0.5, 0.3, 0.01, 0.4])
    end
end

exportgraphics(gcf,'D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Figure\fig3.5chl_az.jpg',...
    'Resolution',600);

%% Fig 3.3.2 Chl, Si, Fe, NH4 in January for 4 experiments and 3 regions 
regionlist = {'saz','pfz','az'};
regionname = {'SAZ','PFZ',' AZ'};
plotnamelist = {'Control','Exp1: D Fe-','Exp5: PI+','Exp6: D Fe- PI+'};
xlim1 = [0,0.5;0,1;0,1.5];
xtick1 = [0.1;0.2;0.3];
xlim2 = [0,5;0,20;10,60];
xtick2 = [1;4;10];
figure(2)
set(gcf,'Position',[100 50 1000 1100])
for i = 1:4
    for j = 1:3
        chl = eval(['TCHL_exp' num2str(i-1) '_' regionlist{j} '(:,7)']);
        si = eval(['Si_exp' num2str(i-1) '_' regionlist{j} '(:,7)']);
        fe = eval(['Fer_exp' num2str(i-1) '_' regionlist{j} '(:,7)'])*1000;
        nh4 = eval(['NH4_exp' num2str(i-1) '_' regionlist{j} '(:,7)']);
        mld = eval(['MLD_' regionlist{j} '(7)']);

        subplot(3,4,4*(j-1)+i)
        set(gca,'Color','none','XTick',[],'YTick',[])
        ax1 = axes('Position',get(gca,'Position'),'NextPlot','add');
        plot(si,depth,'b','LineWidth',1.5,'Parent',ax1);hold on;
        plot([xlim2(j,1) xlim2(j,2)],[mld mld],'Color',[.3 .3 .3],'LineStyle','--')
        set(ax1,'YDir','reverse')
        xlabel(ax1,'Si (\muM)','FontSize',10)
        if i == 1
            ylabel(ax1,'Depth (m)','FontSize',10)
        end
        xlim([xlim2(j,1) xlim2(j,2)])
        xticks(xlim2(j,1):xtick2(j):xlim2(j,2))
        ylim([-inf 300])
        if i ~= 1
            yticklabels([])
        end

        ax1.XAxisLocation = 'top';

        ax2 = axes('Position',get(gca,'Position'),'Color','none','NextPlot','add');
        plot(chl,depth,'Color',[0 .5 0],'LineWidth',1.5,'Parent',ax2);hold on;
        plot([-2,-1],[100,101],'b','LineWidth',1.5,'Parent',ax2)
        plot(fe,depth,'r','LineWidth',1.5,'Parent',ax2)
        plot(nh4,depth,'Color',[.8 .8 0],'LineWidth',1.5,'Parent',ax2)
        set(ax2,'YDir','reverse')
        xlabel(ax2,'Chl (\mug/L), NH4 (\muM), Fe (nM)')
        ylim([-inf 300])
        xlim([xlim1(j,1) xlim1(j,2)])
        xticks(xlim1(j,1):xtick1(j):xlim1(j,2))
        yticklabels([])
        ax2.XAxisLocation = 'bottom';

        set(ax2,'Position',get(ax1,'Position'))

        if i == 4 & j == 1
            legend(ax2,{'Chl','Si','Fe','NH4'},'Color','w','Location','southeast','FontSize',10);
        end

        if j == 1
            t = title(plotnamelist{i},"FontSize",10);
            set(t,'Position',[0.25,-45,0])
        end

        if i == 1
            text(ax2,-0.35*xlim1(j,2),170,regionname{j},"FontSize",10,"FontWeight","bold","Rotation",90)
        end

        text(ax2,0,-40,[char('a'+4*j+i-5) ')'],'FontSize',10,'FontWeight','bold')

        box on;
    end
end

exportgraphics(gcf,'D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Figure\fig3.6var_exp_reg.jpg',...
    'Resolution',600);

%% Fig 4.1 C:Chl ratio in January for 4 experiments and 3 regions 
collist = {'k','g','y','m'};
regionname = {'a) SAZ','b) PFZ','c) AZ'};
figure(3)
set(gcf,'Position',[100 100 1100 400])
for i = 1:length(regionlist)
    subplot(1,4,i)
    hold on;
    for j = 1:length(explist)
        data = eval(['rCChl_exp' num2str(j-1) '_' regionlist{i} '(:,7)']);
        plot(data,depth,collist{j},"LineWidth",1.5)
    end
    mld = eval(['MLD_' regionlist{i} '(7)']);
    plot([0 200],[mld mld],'Color',[.3 .3 .3],'LineStyle','--')
    if i == 1
        lgd=legend(plotnamelist,'FontSize',10);
        set(lgd,'Position',[0.8,0.5,0,0])
        ylabel('Depth (m)','FontSize',10)
    else
        yticklabels([])
    end
    xlim([0 200])
    ylim([-inf 300])
    title(regionname{i},'FontSize',10)
    xlabel('C:Chl (g:g)','FontSize',10)
    set(gca,'XAxisLocation','top','YDir','reverse')
    box on;
end

exportgraphics(gcf,'D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Figure\fig4.1rCChl.jpg',...
    'Resolution',600);



function output=shyear(lastyr,nextyr)
% combine two years datasets and keep July to December of last year and
% January to June of next year
output = lastyr;
output(:,:,:,1:6) = lastyr(:,:,:,7:12);
output(:,:,:,7:12) = nextyr(:,:,:,1:6);
end

function output=shyear2(lastyr,nextyr)
% combine two years datasets and keep July to December of last year and
% January to June of next year
output = lastyr;
output(:,:,1:6) = lastyr(:,:,7:12);
output(:,:,7:12) = nextyr(:,:,1:6);
end