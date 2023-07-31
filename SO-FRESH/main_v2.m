%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to read and do the validation of SOFRESH SSS products using
% Didac Costa's in-situ measurements
% 0. load fronts
% 1. read bec products
% 2. read in-situs
% 3. regrid bec products
% 4. collocalize in-situ to satellite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. loads Fronts Matlab file
load fronts.mat

%% 1. Read BEC Products

pathdir = '/Users/mumbert/Library/Mobile Documents/com~apple~CloudDocs/Analysis/Validation/Validation_SOFRESH/SOFRESH_data/';

fid = fopen([pathdir,'list_l3.txt']);

C   = textscan(fid,'%s');
fclose(fid);
filename = cell2mat(C{1}(1));
bec.lat(:,:)       = ncread(sprintf('%s%s',pathdir,filename),'lat');
bec.lon(:,:)       = ncread(sprintf('%s%s',pathdir,filename),'lon');

for it=1:numel(C{1})
    filename = cell2mat(C{1}(it));
    filedate = datenum(filename(end-35:end-28),'yyyymmdd');
    bec.time(it)       = filedate;
    bec.sss(:,:,it)    = ncread(sprintf('%s%s',pathdir,filename),'sss');
end

%% 2. Read in-situ

load bwr_2011.dat

lat_2011=bwr_2011(:,5);
lon_2011 = bwr_2011(:,6);
temp_2011= bwr_2011(:,7);
sal_2011 = bwr_2011(:,8);
year_2011 = bwr_2011(:,3);
day_2011 = bwr_2011(:,1);
month_2011 = bwr_2011(:,2);
time_2011=bwr_2011(:,4);
hour_2011=floor(time_2011);
minute_2011=(time_2011 - hour_2011)*100;
seconds_2011=zeros(size(year_2011));

data_2011=datenum(datetime(year_2011, month_2011, day_2011, hour_2011, floor(minute_2011), seconds_2011));

aa = find(~isnan(sal_2011) & lat_2011<=-30); % Only points southern than 40 south

coloc_sss_nonan = sal_2011(aa);
coloc_lon_nonan = lon_2011(aa);
coloc_lat_nonan = lat_2011(aa);
coloc_time_nonan = data_2011(aa);

%% some figures

figure 
set(gcf, 'position',[149 108 1000 1000], 'color', 'w');
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(bec.lon,bec.lat,bec.sss(:,:,1)); shading flat;
hold on
%situ
a=1:15:length(coloc_lon_nonan);
m_scatter(coloc_lon_nonan(a),coloc_lat_nonan(a),80,coloc_sss_nonan(a),'filled','MarkerEdgeColor','k');
%FRONTS
m_plot(saf(:,1), saf(:,2),'LineStyle','-', 'color',rgb('Black'), 'LineWidth', 1) 
m_plot(stf(:,1), stf(:,2),'LineStyle', '-.','color',rgb('Black'), 'LineWidth', 1) 

c=colorbar;
cmocean('haline')
c.Label.String = ('PSU ');
caxis([32.5 35.5])
m_coast('patch',[.7 .7 .7],'edgecolor','k');
set(gca,'fontsize',18) 
m_grid('xtick',12,'tickdir','out','linest','-','linewidth',3,'fontsize',10,'xaxisloc','top','yaxisloc','right','fontsize',18);
legend('2011','SMOS SSS','in-situ SSS','Subantarctic Front', 'Subtropical Front','Location','westoutside')
title('')

%% 3. REGRID regular of BEC

cellsize = 0.25;
for i = 1:length(bec.time)
    i
    [regrid_sat, refvec3] = geoloc2grid(double(bec.lat),double(bec.lon), double(bec.sss(:,:,i)), cellsize);
    bec_reg(:,:,i) = regrid_sat;
end
 
llat2 = -90:cellsize:81.9; llon2 = -179.9:cellsize:179.9;
[xx2,yy2]=meshgrid(llon2,llat2);
% catds_date = sat_sss.smos_catds_cpdc_L3.date;
 
save bec_reg yy2 xx2 %catds_regrided.mat catds_reg xx2 yy2 catds_date -v7.3
%load catds_regrided.mat

%sat = bec_reg;
%lon_jpl = xx2;
%lat_jpl = yy2;
%jpl_t = catds_date;


%% 4. Colocation Satellite - BWR 2011
clear BWR_INT ANOM_bec BWR_INTERP RMS_int BWR_RMS SAT_bec
BWR_INT = zeros(688,1440,198);
ANOM_bec = zeros(688,1440,198);
SAT_bec = zeros(688,1440,198);
RMS_int = zeros(688,1440,198);

for id=1:1 %numel(C{1})
    id   
    t0=bec.time(id);t1=t0 +1;
    ind = find(coloc_time_nonan>=t0 & coloc_time_nonan<=t1);
    if ~isempty(ind)
         disp('full')
         bwr_lon = double(coloc_lon_nonan(ind));bwr_lat = double(coloc_lat_nonan(ind));bwr_sss = double(coloc_sss_nonan(ind));
        %TSG_INTERP(:,:,id) = ffgridxvyv(tsg_lon,tsg_lat,tsg_sss,smos_lon,smos_lat);
        %[TSG_INTERP(:,:,id),TSG_RMS(:,:,id), xxvec, yyvec,ngrid] = ffgridrms(tsg_lon,tsg_lat,tsg_sss, 0.25, 0.25, min(smos_lon), min(smos_lat), max(smos_lon), max(smos_lat));
        [BWR_INTERP,BWR_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(bwr_lon),squeeze(bwr_lat),squeeze(bwr_sss)',0.25,0.25,min(xx2(:,1)),min(yy2(1,:)),max(xx2(:,1)),max(yy2(1,:)));
        BWR_INT(:,:,id) = BWR_INTERP;
        SAT_MAP = squeeze(bec_reg(:,:,id));
        ANOM_bec(:,:,id) = SAT_MAP - BWR_INTERP;
        SAT_bec(:,:,id) = SAT_MAP - BWR_INTERP + BWR_INTERP;
        RMS_int(:,:,id) = BWR_RMS;
        non = find(~isnan(BWR_INTERP));
        clear BWR_INTERP BWR_RMS

    else
        disp('empty')
        BWR_INT(:,:,id) = nan;
        ANOM_bec(:,:,id) = nan;
        RMS_int(:,:,id) = nan;
        SAT_bec(:,:,id) = nan;
    end
end
% 
BWR_MEAN = nanmean(BWR_INT,3);
ANOM_MEAN = nanmean(ANOM_bec,3); 
RMS_MEAN = nanmean(RMS_int,3); 
SAT_MEAN = nanmean(SAT_bec,3);

plot_SOFRESH_drake(xx2,yy2,BWR_MEAN,'PSU');caxis([32.5 34.4]);
% %plot_VG(lon_jpl,lat_jpl,ANOM_MEAN','')
% %plot_VG(lon_jpl,lat_jpl,RMS_MEAN','')
% plot_VG_drake(lon_jpl,lat_jpl,SAT_MEAN','');caxis([32.5 34.4]);
% 

% figure 
% set(gcf, 'position',[149 108 1000 1000], 'color', 'w');
% m_proj('stereographic','lat',-90,'long',0,'radius',60);
% m_pcolor(xx2,yy2,bec_reg(:,:,1)); shading flat;
% hold on
% %situ
% a=1:15:length(coloc_lon_nonan);
% m_scatter(coloc_lon_nonan(a),coloc_lat_nonan(a),80,coloc_sss_nonan(a),'filled','MarkerEdgeColor','k');
% %FRONTS
% m_plot(saf(:,1), saf(:,2),'LineStyle','-', 'color',rgb('Black'), 'LineWidth', 1) 
% m_plot(stf(:,1), stf(:,2),'LineStyle', '-.','color',rgb('Black'), 'LineWidth', 1) 
% 
% c=colorbar;
% cmocean('haline')
% c.Label.String = ('PSU ');
% caxis([32.5 35.5])
% m_coast('patch',[.7 .7 .7],'edgecolor','k');
% set(gca,'fontsize',18) 
% m_grid('xtick',12,'tickdir','out','linest','-','linewidth',3,'fontsize',10,'xaxisloc','top','yaxisloc','right','fontsize',18);
% legend('2011','SMOS SSS','in-situ SSS','Subantarctic Front', 'Subtropical Front','Location','westoutside')
% title('')
% % 



