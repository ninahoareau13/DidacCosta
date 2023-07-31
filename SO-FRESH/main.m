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

% 3. 

%% 4. Colocation Satellite - BWR 2011
clear BWR_INT ANOM_bec BWR_INTERP RMS_int BWR_RMS SAT_bec
BWR_INT = zeros(720,720,198);
ANOM_bec = zeros(720,720,198);
SAT_bec = zeros(720,720,198);
RMS_int = zeros(720,720,198);

return
for id=1:1 %numel(C{1})
    id   
    t0=bec.time(id);t1=t0 +1;
    ind = find(coloc_time_nonan>=t0 & coloc_time_nonan<=t1);
    if ~isempty(ind)
         disp('full')
         bwr_lon = double(coloc_lon_nonan(ind));bwr_lat = double(coloc_lat_nonan(ind));bwr_sss = double(coloc_sss_nonan(ind));
        %TSG_INTERP(:,:,id) = ffgridxvyv(tsg_lon,tsg_lat,tsg_sss,smos_lon,smos_lat);
        %[TSG_INTERP(:,:,id),TSG_RMS(:,:,id), xxvec, yyvec,ngrid] = ffgridrms(tsg_lon,tsg_lat,tsg_sss, 0.25, 0.25, min(smos_lon), min(smos_lat), max(smos_lon), max(smos_lat));
        [BWR_INTERP,BWR_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(bwr_lon),squeeze(bwr_lat),squeeze(bwr_sss)',0.25,0.25,min(bec.lon(:,1)),min(bec.lat(1,:)),max(bec.lon(:,1)),max(bec.lat(1,:)));
        BWR_INT(:,:,id) = BWR_INTERP;
        SAT_MAP = squeeze(bec.sss(:,:,id))';
        ANOM_bec(:,:,id) = SAT_MAP - BWR_INTERP;
        SAT_bec(:,:,id) = SAT_MAP - BWR_INTERP + BWR_INTERP;
        RMS_int(:,:,id) = BWR_RMS;
        non = find(~isnan(BWR_INTERP));
        clear BWR_INTERP BWR_RMS
    end
end



%     else
%         disp('empty')
%         TSG_INT(:,:,id) = nan;
%         ANOM_jpl(:,:,id) = nan;
%         RMS_int(:,:,id) = nan;
%         SAT_jpl(:,:,id) = nan;
%     end
% end
% 
% TSG_MEAN = nanmean(TSG_INT,3);ANOM_MEAN = nanmean(ANOM_jpl,3); RMS_MEAN = nanmean(RMS_int,3); SAT_MEAN = nanmean(SAT_jpl,3);
% 
% plot_VG_drake(lon_jpl,lat_jpl,TSG_MEAN','');caxis([32.5 34.4]);
% %plot_VG(lon_jpl,lat_jpl,ANOM_MEAN','')
% %plot_VG(lon_jpl,lat_jpl,RMS_MEAN','')
% plot_VG_drake(lon_jpl,lat_jpl,SAT_MEAN','');caxis([32.5 34.4]);
% 


