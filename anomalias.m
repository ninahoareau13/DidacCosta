clear all

% ANOMALIES

% Obejtivo: Optimizar programa actual
% Climatology Salinity WOA 2013
sal_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "s_an"); 
sal_woa = sal_woa(:,:,1); %coger superficie
lon_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "lon");
lat_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "lat"); 
[matriz_lat_woa,matriz_lon_woa]=meshgrid(lat_woa,lon_woa);
jpl_t=ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "time");
jpl=sal_woa;

% load Polar Fronts data
load('C:\Users\Usuario\Desktop\TFM\Datos complementarios\fronts.mat');


%% 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sea Ice concentration 2020
ice_conc_2020 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_20202021_NOVFEB_mean.nc", "ice_conc"); 
ice_lat_2020= ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_20202021_NOVFEB_mean.nc", "lat");
ice_lon_2020 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_20202021_NOVFEB_mean.nc", "lon");

% load VG data 2020
load('C:\Users\Usuario\Desktop\TFM\Datos regatas\JMSE_data_filtered_VG.mat')
% 
% ff_2020=find(Lat_mat<=-40);
% 
% lat_2020=Lat_mat(ff_2020);
% lon_2020=Lon_mat(ff_2020);
% sal_2020=Sal_mat(ff_2020);
% temp_2020 = Temp_mat(ff_2020);
% tt=Data_mat(ff_2020);

%% SALINIDAD 
%%COLOCATION SATELITE  - Vendee Globe

aa = find(~isnan(Sal_mat) & Lat_mat<=-40); % Only points southern than 40 south

coloc_sss_nonan = Sal_mat(aa);
coloc_lon_nonan = Lon_mat(aa);
coloc_lat_nonan = Lat_mat(aa);
coloc_time_nonan = Data_mat(aa);

%%COLOCATION WOA

% clear TSG_INT ANOM_jpl TSG_INTERP RMS_int TSG_RMS SAT_jpl
% TSG_INT = zeros(720,1440,length(jpl_t)); % SMAP
% ANOM_jpl = zeros(720,1440,length(jpl_t));
% SAT_jpl = zeros(720,1440,length(jpl_t));
% RMS_int = zeros(720,1440,length(jpl_t));

% TSG_INT = zeros(720,1440,1); % 1 porque WOA solo tiene una data
% ANOM_jpl = zeros(720,1440,1);
% SAT_jpl = zeros(720,1440,1);
% RMS_int = zeros(720,1440,1);

TSG_INT = zeros(720,1440,1); % 1  WOA date
ANOM_jpl = zeros(720,1440,1);
SAT_jpl = zeros(720,1440,1);
RMS_int = zeros(720,1440,1);

% TSG_INT = zeros(3600,7200,length(jpl_t)); % BEC L4
% ANOM_jpl = zeros(3600,7200,length(jpl_t));
% RMS_int = zeros(3600,7200,length(jpl_t));


tsg_lon = double(coloc_lon_nonan);tsg_lat = double(coloc_lat_nonan);tsg_sss = double(coloc_sss_nonan);tsg_data=double(coloc_time_nonan);
%TSG_INTERP(:,:,id) = ffgridxvyv(tsg_lon,tsg_lat,tsg_sss,smos_lon,smos_lat);
 %[TSG_INTERP(:,:,id),TSG_RMS(:,:,id), xxvec, yyvec,ngrid] = ffgridrms(tsg_lon,tsg_lat,tsg_sss, 0.25, 0.25, min(smos_lon), min(smos_lat), max(smos_lon), max(smos_lat));
[TSG_INTERP,TSG_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_sss)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
[TSG_INTERP_data,TSG_RMS_data,xxvec, yyvec_data,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_data)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
%donde sss meter tiempo (tsg:tiempo)
% [TSG_INTERP,TSG_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_sss),0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
        
TSG_INT(:,:,1) = TSG_INTERP;   % 1 - una data
SAT_MAP = squeeze(jpl(:,:,1))';
% ANOM_jpl(:,:,1) = SAT_MAP - TSG_INTERP;
ANOM_jpl(:,:,1) = -(SAT_MAP - TSG_INTERP);
SAT_jpl(:,:,1) = SAT_MAP - TSG_INTERP + TSG_INTERP;
RMS_int(:,:,1) = TSG_RMS;
non = find(~isnan(TSG_INTERP));
%clear TSG_INTERP TSG_RMS

[matriz_lon_woa,matriz_lat_woa]=meshgrid(lon_woa,lat_woa);

% FIGURE...............................................................
figure 
set(gcf,'color','w');
% m_proj('ortho','lat',30,'long',-20');
% m_grid('linest','-','xticklabels',[],'yticklabels',[]);
m_proj('stereographic','lat',-90,'long',0,'radius',50);

%FRONTS
hold on
m_plot(saf(:,1), saf(:,2),'LineStyle','-', 'color',rgb("Black"), 'LineWidth', 1) 
m_plot(stf(:,1), stf(:,2),'LineStyle', '-.','color',rgb("Black"), 'LineWidth', 1) 
%.................................
%ICE CONC........................
hold on; v=[15,15]; m_contour(ice_lon_2020, ice_lat_2020, ice_conc_2020, v, '-', 'color',rgb('Blue'), 'Linewidth', 2)
%................................

% ff=find(~isnan(ANOM_MEAN));
ff=find(~isnan(ANOM_jpl));
latt_woa=matriz_lat_woa(ff);
lonn_woa=matriz_lon_woa(ff);
% anomalies=ANOM_MEAN(ff);
anomalies=ANOM_jpl(ff);
a=1:5:length(lonn_woa);
m_scatter(lonn_woa(a),latt_woa(a),80,anomalies(a),'filled','MarkerEdgeColor','k');
c=colorbar;
cmocean('balance')
c.Label.String = ('Salinity Anomalies');
caxis([-1.5 1.5])
m_coast('patch',[.7 .7 .7],'edgecolor','k');
set(gca,'fontsize',18) 
m_grid('xtick',12,'tickdir','out','linest','-','linewidth',3,'fontsize',10,'xaxisloc','top','yaxisloc','right','fontsize',18);
%title('Salinity')
legend('', 'Subantarctic Front', 'Subtropical Front', 'Ice Conc. 15%', 'Location','westoutside')

% Serie temporal
d=find(~isnan(ANOM_jpl));
data_interp=datetime(TSG_INTERP_data, 'ConvertFrom', 'datenum');

figure
set(gcf,'color','w');
scatter(data_interp(d),ANOM_jpl(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
c=colorbar;
cmocean('balance')
% colormap(redbluecmap)
c.Label.String = ('Anomalies (PSU)');
caxis([-1.5 1.5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Anomalies (PSU)')
title('Salinity Anomalies')
subtitle('Southern Ocean - 2020', 'FontSize', 14);
grid on
%xline(data, '--r')

figure
set(gcf,'color','w');
scatter(data_interp(d),TSG_INTERP(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
c=colorbar;
cmocean('balance')
c.Label.String = ('Anomalies (PSU)');
caxis([-1.5 1.5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Salinity (PSU)')
ylim([32 36])
title('Salinity Anomalies')
subtitle('Southern Ocean - 2020', 'FontSize', 14);
grid on
%xline

%% TEMPERATURE

clear all

% Sea Ice concentration 2020
ice_conc_2020 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_20202021_NOVFEB_mean.nc", "ice_conc"); 
ice_lat_2020= ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_20202021_NOVFEB_mean.nc", "lat");
ice_lon_2020 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_20202021_NOVFEB_mean.nc", "lon");

% load VG data 2020
load('C:\Users\Usuario\Desktop\TFM\Datos regatas\JMSE_data_filtered_VG.mat')

% ANOMALIES

% Climatology WOA 2013
temp_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "t_an");
temp_woa = temp_woa(:,:,1);
lon_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "lon");
lat_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "lat");
[matriz_lat_woa,matriz_lon_woa]=meshgrid(lat_woa,lon_woa);
jpl_t=ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "time");
jpl=temp_woa;

% load Polar Fronts data
load('C:\Users\Usuario\Desktop\TFM\Datos complementarios\fronts.mat');
%.................................................................................

aa = find(~isnan(Temp_mat) & Lat_mat<=-40); % Only points southern than 40 south

coloc_sst_nonan = Temp_mat(aa);
coloc_lon_nonan = Lon_mat(aa);
coloc_lat_nonan = Lat_mat(aa);
coloc_time_nonan = Data_mat(aa);

TSG_INT = zeros(720,1440,1); % 1  WOA date
ANOM_jpl = zeros(720,1440,1);
SAT_jpl = zeros(720,1440,1);
RMS_int = zeros(720,1440,1);

tsg_lon = double(coloc_lon_nonan);tsg_lat = double(coloc_lat_nonan);tsg_sst = double(coloc_sst_nonan);tsg_data=double(coloc_time_nonan);
[TSG_INTERP,TSG_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_sst)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
[TSG_INTERP_data,TSG_RMS_data,xxvec, yyvec_data,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_data)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
 
TSG_INT(:,:,1) = TSG_INTERP;   % 1 porque solo hay una data
SAT_MAP = squeeze(jpl(:,:,1))';
ANOM_jpl(:,:,1) = -(SAT_MAP - TSG_INTERP);
SAT_jpl(:,:,1) = SAT_MAP - TSG_INTERP + TSG_INTERP;
RMS_int(:,:,1) = TSG_RMS;
non = find(~isnan(TSG_INTERP));
%clear TSG_INTERP TSG_RMS

[matriz_lon_woa,matriz_lat_woa]=meshgrid(lon_woa,lat_woa);

% FIGURE...............................................................
figure 
set(gcf,'color','w');
m_proj('stereographic','lat',-90,'long',0,'radius',50);

%FRONTS
hold on
m_plot(saf(:,1), saf(:,2),'LineStyle','-', 'color',rgb("Black"), 'LineWidth', 1) 
m_plot(stf(:,1), stf(:,2),'LineStyle', '-.','color',rgb("Black"), 'LineWidth', 1) 
%.................................
%ICE CONC........................
hold on; v=[15,15]; m_contour(ice_lon_2020, ice_lat_2020, ice_conc_2020, v, '-', 'color',rgb('Blue'), 'Linewidth', 2)
%................................

% ff=find(~isnan(ANOM_MEAN));
ff=find(~isnan(ANOM_jpl));
latt_woa=matriz_lat_woa(ff);
lonn_woa=matriz_lon_woa(ff);
% anomalies=ANOM_MEAN(ff);
anomalies=ANOM_jpl(ff);
a=1:5:length(lonn_woa);
m_scatter(lonn_woa(a),latt_woa(a),80,anomalies(a),'filled','MarkerEdgeColor','k');
c=colorbar;
cmocean('balance')
c.Label.String = ('Temperature Anomalies (ºC)');
caxis([-5 5])
m_coast('patch',[.7 .7 .7],'edgecolor','k');
set(gca,'fontsize',18) 
m_grid('xtick',12,'tickdir','out','linest','-','linewidth',3,'fontsize',10,'xaxisloc','top','yaxisloc','right','fontsize',18);
%title('Salinity')
legend('', 'Subantarctic Front', 'Subtropical Front', 'Ice Conc. 15%', 'Location','westoutside')

% Serie temporal
d=find(~isnan(ANOM_jpl));
data_interp=datetime(TSG_INTERP_data, 'ConvertFrom', 'datenum');

figure
set(gcf,'color','w');
scatter(data_interp(d),ANOM_jpl(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
c=colorbar;
cmocean('balance')
caxis([-5 5])
% colormap(redbluecmap)
c.Label.String = ('Anomalies (ºC)');
%caxis([-1.5 1])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Anomalies (ºC)')
title('Temperature Anomalies')
subtitle('Southern Ocean - 2020', 'FontSize', 14);
grid on
%xline(data, '--r')

figure
set(gcf,'color','w');
scatter(data_interp(d),TSG_INTERP(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
ylim([4 22])
c=colorbar;
cmocean('balance')
c.Label.String = ('Anomalies (ºC)');
caxis([-5 5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Temperature (ºC)')
title('Temperature Anomalies')
subtitle('Southern Ocean - 2020', 'FontSize', 14);
grid on


%% 2011

% Climatology Salinity WOA 2013
sal_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "s_an"); 
sal_woa = sal_woa(:,:,1); %coger superficie
lon_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "lon");
lat_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "lat"); 
[matriz_lat_woa,matriz_lon_woa]=meshgrid(lat_woa,lon_woa);
jpl_t=ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "time");
jpl=sal_woa;

% load Polar Fronts data
load('C:\Users\Usuario\Desktop\TFM\Datos complementarios\fronts.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sea Ice concentration 2011
ice_conc_2011 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2011_JANAPR_mean.nc", "ice_conc"); 
ice_lat_2011 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2011_JANAPR_mean.nc", "lat");
ice_lon_2011 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2011_JANAPR_mean.nc", "lon");

% load BWR data - 2011
load('C:\Users\Usuario\Desktop\TFM\Datos regatas\bwr_2011.dat')

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


%%SALINIDAD 
%%COLOCATION SATELITE  - Vendee Globe

aa = find(~isnan(sal_2011) & lat_2011<=-40); % Only points southern than 40 south

coloc_sss_nonan = sal_2011(aa);
coloc_lon_nonan = lon_2011(aa);
coloc_lat_nonan = lat_2011(aa);
coloc_time_nonan = data_2011(aa);

%%COLOCATION WOA

TSG_INT = zeros(720,1440,1); % 1  WOA date
ANOM_jpl = zeros(720,1440,1);
SAT_jpl = zeros(720,1440,1);
RMS_int = zeros(720,1440,1);

% TSG_INT = zeros(3600,7200,length(jpl_t)); % BEC L4
% ANOM_jpl = zeros(3600,7200,length(jpl_t));
% RMS_int = zeros(3600,7200,length(jpl_t));


tsg_lon = double(coloc_lon_nonan);tsg_lat = double(coloc_lat_nonan);tsg_sss = double(coloc_sss_nonan);tsg_data=double(coloc_time_nonan);
%TSG_INTERP(:,:,id) = ffgridxvyv(tsg_lon,tsg_lat,tsg_sss,smos_lon,smos_lat);
 %[TSG_INTERP(:,:,id),TSG_RMS(:,:,id), xxvec, yyvec,ngrid] = ffgridrms(tsg_lon,tsg_lat,tsg_sss, 0.25, 0.25, min(smos_lon), min(smos_lat), max(smos_lon), max(smos_lat));
[TSG_INTERP,TSG_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_sss)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
[TSG_INTERP_data,TSG_RMS_data,xxvec, yyvec_data,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_data)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
%donde sss meter tiempo (tsg:tiempo)
% [TSG_INTERP,TSG_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_sss),0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
        
TSG_INT(:,:,1) = TSG_INTERP;   % 1 porque solo hay una data
SAT_MAP = squeeze(jpl(:,:,1))';
% ANOM_jpl(:,:,1) = SAT_MAP - TSG_INTERP;
ANOM_jpl(:,:,1) = -(SAT_MAP - TSG_INTERP);
SAT_jpl(:,:,1) = SAT_MAP - TSG_INTERP + TSG_INTERP;
RMS_int(:,:,1) = TSG_RMS;
non = find(~isnan(TSG_INTERP));
%clear TSG_INTERP TSG_RMS

%TSG_MEAN = nanmean(TSG_INT,3);ANOM_MEAN = nanmean(ANOM_jpl,3); RMS_MEAN = nanmean(RMS_int,3); SAT_MEAN = nanmean(SAT_jpl,3);

% plot_VG(lon_jpl,lat_jpl,TSG_MEAN,'')
% plot_VG(lon_jpl,lat_jpl,ANOM_MEAN,'')
% plot_VG(lon_jpl,lat_jpl,RMS_MEAN,'')

[matriz_lon_woa,matriz_lat_woa]=meshgrid(lon_woa,lat_woa);

% FIGURE...............................................................
figure 
set(gcf,'color','w');
% m_proj('ortho','lat',30,'long',-20');
% m_grid('linest','-','xticklabels',[],'yticklabels',[]);
m_proj('stereographic','lat',-90,'long',0,'radius',50);

%FRONTS...........................
hold on
m_plot(saf(:,1), saf(:,2),'LineStyle','-', 'color',rgb("Black"), 'LineWidth', 1) 
m_plot(stf(:,1), stf(:,2),'LineStyle', '-.','color',rgb("Black"), 'LineWidth', 1) 
%.................................
%ICE CONC........................
hold on; v=[15,15]; m_contour(ice_lon_2011, ice_lat_2011, ice_conc_2011, v, '-', 'color',rgb('Blue'), 'Linewidth', 2)
%................................

% ff=find(~isnan(ANOM_MEAN));
ff=find(~isnan(ANOM_jpl));
latt_woa=matriz_lat_woa(ff);
lonn_woa=matriz_lon_woa(ff);
% anomalies=ANOM_MEAN(ff);
anomalies=ANOM_jpl(ff);
a=1:5:length(lonn_woa);
m_scatter(lonn_woa(a),latt_woa(a),80,anomalies(a),'filled','MarkerEdgeColor','k');
%m_pcolor(matriz_lon_woa,matriz_lat_woa,40,ANOM_MEAN,'filled');
%m_pcolor(matriz_lon_woa,matriz_lat_woa,ANOM_MEAN,'MarkerSize',40)
c=colorbar;
cmocean('balance')
c.Label.String = ('Salinity Anomalies (PSU)');
caxis([-1.5 1.5])
m_coast('patch',[.7 .7 .7],'edgecolor','k');
set(gca,'fontsize',18) 
m_grid('xtick',12,'tickdir','out','linest','-','linewidth',3,'fontsize',10,'xaxisloc','top','yaxisloc','right','fontsize',18);
%title('Salinity')
legend('', 'Subantarctic Front', 'Subtropical Front', 'Ice Conc. 15%', 'Location','westoutside')

% Serie temporal
d=find(~isnan(ANOM_jpl));
data_interp=datetime(TSG_INTERP_data, 'ConvertFrom', 'datenum');

figure
set(gcf,'color','w');
scatter(data_interp(d),ANOM_jpl(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
c=colorbar;
cmocean('balance')
% colormap(redbluecmap)
c.Label.String = ('Anomalies (PSU)');
caxis([-1.5 1.5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Anomalies (PSU)')
title('Salinity Anomalies')
subtitle('Southern Ocean - 2011', 'FontSize', 14);
grid on
%xline(data, '--r')

figure
set(gcf,'color','w');
scatter(data_interp(d),TSG_INTERP(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
ylim([32 36])
c=colorbar;
cmocean('balance')
c.Label.String = ('Anomalies (PSU)');
caxis([-1.5 1.5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Salinity (PSU)')
title('Salinity Anomalies')
subtitle('Southern Ocean - 2011', 'FontSize', 14);
grid on

%% TEMPERATURE

% Climatology Salinity WOA 2013
temp_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "t_an");
temp_woa = temp_woa(:,:,1);
lon_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "lon");
lat_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "lat");
[matriz_lat_woa,matriz_lon_woa]=meshgrid(lat_woa,lon_woa);
jpl_t=ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "time");
jpl=temp_woa;

% load Polar Fronts data
load('C:\Users\Usuario\Desktop\TFM\Datos complementarios\fronts.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sea Ice concentration 2011
ice_conc_2011 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2011_JANAPR_mean.nc", "ice_conc"); 
ice_lat_2011 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2011_JANAPR_mean.nc", "lat");
ice_lon_2011 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2011_JANAPR_mean.nc", "lon");

% load BWR data - 2011
load('C:\Users\Usuario\Desktop\TFM\Datos regatas\bwr_2011.dat')

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

%%COLOCATION SATELITE  - Vendee Globe

aa = find(~isnan(temp_2011) & lat_2011<=-40); % Only points southern than 40 south

coloc_sst_nonan = temp_2011(aa);
coloc_lon_nonan = lon_2011(aa);
coloc_lat_nonan = lat_2011(aa);
coloc_time_nonan = data_2011(aa);

TSG_INT = zeros(720,1440,1); % 1  WOA date
ANOM_jpl = zeros(720,1440,1);
SAT_jpl = zeros(720,1440,1);
RMS_int = zeros(720,1440,1);

tsg_lon = double(coloc_lon_nonan);tsg_lat = double(coloc_lat_nonan);tsg_sst = double(coloc_sst_nonan);tsg_data=double(coloc_time_nonan);
[TSG_INTERP,TSG_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_sst)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
[TSG_INTERP_data,TSG_RMS_data,xxvec, yyvec_data,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_data)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
 
TSG_INT(:,:,1) = TSG_INTERP;   % 1 porque solo hay una data
SAT_MAP = squeeze(jpl(:,:,1))';
ANOM_jpl(:,:,1) = -(SAT_MAP - TSG_INTERP);
SAT_jpl(:,:,1) = SAT_MAP - TSG_INTERP + TSG_INTERP;
RMS_int(:,:,1) = TSG_RMS;
non = find(~isnan(TSG_INTERP));
%clear TSG_INTERP TSG_RMS

[matriz_lon_woa,matriz_lat_woa]=meshgrid(lon_woa,lat_woa);

% FIGURE...............................................................
figure 
set(gcf,'color','w');
m_proj('stereographic','lat',-90,'long',0,'radius',50);

%FRONTS
hold on
m_plot(saf(:,1), saf(:,2),'LineStyle','-', 'color',rgb("Black"), 'LineWidth', 1) 
m_plot(stf(:,1), stf(:,2),'LineStyle', '-.','color',rgb("Black"), 'LineWidth', 1) 
%.................................
%ICE CONC........................
hold on; v=[15,15]; m_contour(ice_lon_2011, ice_lat_2011, ice_conc_2011, v, '-', 'color',rgb('Blue'), 'Linewidth', 2)
%................................

% ff=find(~isnan(ANOM_MEAN));
ff=find(~isnan(ANOM_jpl));
latt_woa=matriz_lat_woa(ff);
lonn_woa=matriz_lon_woa(ff);
% anomalies=ANOM_MEAN(ff);
anomalies=ANOM_jpl(ff);
a=1:5:length(lonn_woa);
m_scatter(lonn_woa(a),latt_woa(a),80,anomalies(a),'filled','MarkerEdgeColor','k');
c=colorbar;
cmocean('balance')
c.Label.String = ('Temperature Anomalies (ºC)');
caxis([-5 5])
m_coast('patch',[.7 .7 .7],'edgecolor','k');
set(gca,'fontsize',18) 
m_grid('xtick',12,'tickdir','out','linest','-','linewidth',3,'fontsize',10,'xaxisloc','top','yaxisloc','right','fontsize',18);
%title('Salinity')
legend('', 'Subantarctic Front', 'Subtropical Front', 'Ice Conc. 15%', 'Location','westoutside')

% Serie temporal
d=find(~isnan(ANOM_jpl));
data_interp=datetime(TSG_INTERP_data, 'ConvertFrom', 'datenum');

figure
set(gcf,'color','w');
scatter(data_interp(d),ANOM_jpl(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
c=colorbar;
cmocean('balance')
%caxis([-5 6])
% colormap(redbluecmap)
c.Label.String = ('Anomalies (ºC)');
caxis([-5 5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Anomalies (ºC)')
title('Temperature Anomalies')
subtitle('Southern Ocean - 2011', 'FontSize', 14);
grid on
%xline(data, '--r')

figure
set(gcf,'color','w');
scatter(data_interp(d),TSG_INTERP(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
ylim([4 22])
c=colorbar;
cmocean('balance')
c.Label.String = ('Anomalies (ºC)');
caxis([-5 5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Temperature (ºC)')
title('Temperature Anomalies')
subtitle('Southern Ocean - 2011', 'FontSize', 14);
grid on
%xline(data, '--r')

%% 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% Climatology Salinity WOA 2013
sal_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "s_an"); 
sal_woa = sal_woa(:,:,1); %coger superficie
lon_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "lon");
lat_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "lat"); 
[matriz_lat_woa,matriz_lon_woa]=meshgrid(lat_woa,lon_woa);
jpl_t=ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_s00_04.nc", "time");
jpl=sal_woa;

% load Polar Fronts data
load('C:\Users\Usuario\Desktop\TFM\Datos complementarios\fronts.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sea Ice concentration 2015
ice_conc_2015 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2015_JANAPR_mean.nc", "ice_conc"); 
ice_lat_2015 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2015_JANAPR_mean.nc", "lat");
ice_lon_2015 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2015_JANAPR_mean.nc", "lon");

% load BWR data 2015
load('C:\Users\Usuario\Desktop\TFM\Datos regatas\Collocated_BWR_2014_2015_GPS_only_SBE_frequency_QC_manual_final_3.mat')

lat_2015=data_tsgqc.LATX;
lon_2015 = data_tsgqc.LONX;
sal_2015=data_tsgqc.SSPS;
data_2015=datenum(datetime(data_tsgqc.YEAR, data_tsgqc.MNTH, data_tsgqc.DAYX, data_tsgqc.hh, data_tsgqc.mi, data_tsgqc.ss));

ffsg_2015=find(data_tsgqc.SSPS_QC==1); % Good data temperature 2015
lat_2015=lat_2015(ffsg_2015);
lon_2015=lon_2015(ffsg_2015);
sal_2015=sal_2015(ffsg_2015);
data_2015=data_2015(ffsg_2015);


%%SALINIDAD 
%%COLOCATION SATELITE  - Vendee Globe

aa = find(~isnan(sal_2015) & lat_2015<=-40); % Only points southern than 40 south

coloc_sss_nonan = sal_2015(aa);
coloc_lon_nonan = lon_2015(aa);
coloc_lat_nonan = lat_2015(aa);
coloc_time_nonan = data_2015(aa);

%%COLOCATION WOA

TSG_INT = zeros(720,1440,1); % 1  WOA date
ANOM_jpl = zeros(720,1440,1);
SAT_jpl = zeros(720,1440,1);
RMS_int = zeros(720,1440,1);

% TSG_INT = zeros(3600,7200,length(jpl_t)); % BEC L4
% ANOM_jpl = zeros(3600,7200,length(jpl_t));
% RMS_int = zeros(3600,7200,length(jpl_t));


tsg_lon = double(coloc_lon_nonan);tsg_lat = double(coloc_lat_nonan);tsg_sss = double(coloc_sss_nonan);tsg_data=double(coloc_time_nonan);
%TSG_INTERP(:,:,id) = ffgridxvyv(tsg_lon,tsg_lat,tsg_sss,smos_lon,smos_lat);
 %[TSG_INTERP(:,:,id),TSG_RMS(:,:,id), xxvec, yyvec,ngrid] = ffgridrms(tsg_lon,tsg_lat,tsg_sss, 0.25, 0.25, min(smos_lon), min(smos_lat), max(smos_lon), max(smos_lat));
[TSG_INTERP,TSG_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_sss)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
[TSG_INTERP_data,TSG_RMS_data,xxvec, yyvec_data,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_data)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
%donde sss meter tiempo (tsg:tiempo)
% [TSG_INTERP,TSG_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_sss),0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
        
TSG_INT(:,:,1) = TSG_INTERP;   % 1 porque solo hay una data
SAT_MAP = squeeze(jpl(:,:,1))';
% ANOM_jpl(:,:,1) = SAT_MAP - TSG_INTERP;
ANOM_jpl(:,:,1) = -(SAT_MAP - TSG_INTERP);
SAT_jpl(:,:,1) = SAT_MAP - TSG_INTERP + TSG_INTERP;
RMS_int(:,:,1) = TSG_RMS;
non = find(~isnan(TSG_INTERP));
%clear TSG_INTERP TSG_RMS

%TSG_MEAN = nanmean(TSG_INT,3);ANOM_MEAN = nanmean(ANOM_jpl,3); RMS_MEAN = nanmean(RMS_int,3); SAT_MEAN = nanmean(SAT_jpl,3);

% plot_VG(lon_jpl,lat_jpl,TSG_MEAN,'')
% plot_VG(lon_jpl,lat_jpl,ANOM_MEAN,'')
% plot_VG(lon_jpl,lat_jpl,RMS_MEAN,'')

[matriz_lon_woa,matriz_lat_woa]=meshgrid(lon_woa,lat_woa);

% FIGURE...............................................................
figure 
set(gcf,'color','w');
% m_proj('ortho','lat',30,'long',-20');
% m_grid('linest','-','xticklabels',[],'yticklabels',[]);
m_proj('stereographic','lat',-90,'long',0,'radius',50);

%FRONTS...........................
hold on
m_plot(saf(:,1), saf(:,2),'LineStyle','-', 'color',rgb("Black"), 'LineWidth', 1) 
m_plot(stf(:,1), stf(:,2),'LineStyle', '-.','color',rgb("Black"), 'LineWidth', 1) 
%.................................
%ICE CONC........................
hold on; v=[15,15]; m_contour(ice_lon_2015, ice_lat_2015, ice_conc_2015, v, '-', 'color',rgb('Blue'), 'Linewidth', 2)
%................................

% ff=find(~isnan(ANOM_MEAN));
ff=find(~isnan(ANOM_jpl));
latt_woa=matriz_lat_woa(ff);
lonn_woa=matriz_lon_woa(ff);
% anomalies=ANOM_MEAN(ff);
anomalies=ANOM_jpl(ff);
a=1:5:length(lonn_woa);
m_scatter(lonn_woa(a),latt_woa(a),80,anomalies(a),'filled','MarkerEdgeColor','k');
%m_pcolor(matriz_lon_woa,matriz_lat_woa,40,ANOM_MEAN,'filled');
%m_pcolor(matriz_lon_woa,matriz_lat_woa,ANOM_MEAN,'MarkerSize',40)
c=colorbar;
cmocean('balance')
c.Label.String = ('Salinity Anomalies (PSU)');
caxis([-1.5 1.5])
m_coast('patch',[.7 .7 .7],'edgecolor','k');
set(gca,'fontsize',18) 
m_grid('xtick',12,'tickdir','out','linest','-','linewidth',3,'fontsize',10,'xaxisloc','top','yaxisloc','right','fontsize',18);
%title('Salinity')
legend('', 'Subantarctic Front', 'Subtropical Front', 'Ice Conc. 15%', 'Location','westoutside')

% Serie temporal
d=find(~isnan(ANOM_jpl));
data_interp=datetime(TSG_INTERP_data, 'ConvertFrom', 'datenum');

figure
set(gcf,'color','w');
scatter(data_interp(d),ANOM_jpl(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
c=colorbar;
cmocean('balance')
% colormap(redbluecmap)
c.Label.String = ('Anomalies (PSU)');
caxis([-1.5 1.5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Anomalies (PSU)')
title('Salinity Anomalies')
subtitle('Southern Ocean - 2015', 'FontSize', 14);
grid on
%xline(data, '--r')

figure
set(gcf,'color','w');
scatter(data_interp(d),TSG_INTERP(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
ylim([32 36])
c=colorbar;
cmocean('balance')
c.Label.String = ('Anomalies (PSU)');
caxis([-1.5 1.5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Salinity (PSU)')
title('Salinity Anomalies')
subtitle('Southern Ocean - 2015', 'FontSize', 14);
grid on

%% TEMPERATURE

% Climatology Salinity WOA 2013
temp_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "t_an");
temp_woa = temp_woa(:,:,1);
lon_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "lon");
lat_woa = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "lat");
[matriz_lat_woa,matriz_lon_woa]=meshgrid(lat_woa,lon_woa);
jpl_t=ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\woa13_decav_t00_04v2.nc", "time");
jpl=temp_woa;

% load Polar Fronts data
load('C:\Users\Usuario\Desktop\TFM\Datos complementarios\fronts.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sea Ice concentration 2015
ice_conc_2015 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2015_JANAPR_mean.nc", "ice_conc"); 
ice_lat_2015 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2015_JANAPR_mean.nc", "lat");
ice_lon_2015 = ncread("C:\Users\Usuario\Desktop\TFM\Datos complementarios\Sea Ice\ESACCI-SEAICE-L4-SICONC-AMSR_25_2015_JANAPR_mean.nc", "lon");

% load BWR data 2015
load('C:\Users\Usuario\Desktop\TFM\Datos regatas\Collocated_BWR_2014_2015_GPS_only_SBE_frequency_QC_manual_final_3.mat')

lat_2015=data_tsgqc.LATX;
lon_2015 = data_tsgqc.LONX;
temp_2015=data_tsgqc.SSJT;
data_2015=datenum(datetime(data_tsgqc.YEAR, data_tsgqc.MNTH, data_tsgqc.DAYX, data_tsgqc.hh, data_tsgqc.mi, data_tsgqc.ss));

fftg_2015=find(data_tsgqc.SSJT_QC==1); % Good data temperature 2015
lat_2015=lat_2015(fftg_2015);
lon_2015=lon_2015(fftg_2015);
temp_2015=temp_2015(fftg_2015);
data_2015=data_2015(fftg_2015);

%%COLOCATION SATELITE  - Vendee Globe

aa = find(~isnan(temp_2015) & lat_2015<=-40); % Only points southern than 40 south

coloc_sst_nonan = temp_2015(aa);
coloc_lon_nonan = lon_2015(aa);
coloc_lat_nonan = lat_2015(aa);
coloc_time_nonan = data_2015(aa);

TSG_INT = zeros(720,1440,1); % 1  WOA date
ANOM_jpl = zeros(720,1440,1);
SAT_jpl = zeros(720,1440,1);
RMS_int = zeros(720,1440,1);

tsg_lon = double(coloc_lon_nonan);tsg_lat = double(coloc_lat_nonan);tsg_sst = double(coloc_sst_nonan);tsg_data=double(coloc_time_nonan);
[TSG_INTERP,TSG_RMS,xxvec, yyvec,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_sst)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
[TSG_INTERP_data,TSG_RMS_data,xxvec, yyvec_data,ngrid] = ffgridrms(squeeze(tsg_lon),squeeze(tsg_lat),squeeze(tsg_data)',0.25,0.25,min(matriz_lon_woa(:,1)),min(matriz_lat_woa(1,:)),max(matriz_lon_woa(:,1)),max(matriz_lat_woa(1,:)));
 
TSG_INT(:,:,1) = TSG_INTERP;   % 1 porque solo hay una data
SAT_MAP = squeeze(jpl(:,:,1))';
ANOM_jpl(:,:,1) = -(SAT_MAP - TSG_INTERP);
SAT_jpl(:,:,1) = SAT_MAP - TSG_INTERP + TSG_INTERP;
RMS_int(:,:,1) = TSG_RMS;
non = find(~isnan(TSG_INTERP));
%clear TSG_INTERP TSG_RMS

[matriz_lon_woa,matriz_lat_woa]=meshgrid(lon_woa,lat_woa);

% FIGURE...............................................................
figure 
set(gcf,'color','w');
m_proj('stereographic','lat',-90,'long',0,'radius',50);

%FRONTS
hold on
m_plot(saf(:,1), saf(:,2),'LineStyle','-', 'color',rgb("Black"), 'LineWidth', 1) 
m_plot(stf(:,1), stf(:,2),'LineStyle', '-.','color',rgb("Black"), 'LineWidth', 1) 
%.................................
%ICE CONC........................
hold on; v=[15,15]; m_contour(ice_lon_2015, ice_lat_2015, ice_conc_2015, v, '-', 'color',rgb('Blue'), 'Linewidth', 2)
%................................

% ff=find(~isnan(ANOM_MEAN));
ff=find(~isnan(ANOM_jpl));
latt_woa=matriz_lat_woa(ff);
lonn_woa=matriz_lon_woa(ff);
% anomalies=ANOM_MEAN(ff);
anomalies=ANOM_jpl(ff);
a=1:5:length(lonn_woa);
m_scatter(lonn_woa(a),latt_woa(a),80,anomalies(a),'filled','MarkerEdgeColor','k');
c=colorbar;
cmocean('balance')
c.Label.String = ('Temperature Anomalies (ºC)');
caxis([-5 5])
m_coast('patch',[.7 .7 .7],'edgecolor','k');
set(gca,'fontsize',18) 
m_grid('xtick',12,'tickdir','out','linest','-','linewidth',3,'fontsize',10,'xaxisloc','top','yaxisloc','right','fontsize',18);
%title('Salinity')
legend('', 'Subantarctic Front', 'Subtropical Front', 'Ice Conc. 15%', 'Location','westoutside')

% Serie temporal
d=find(~isnan(ANOM_jpl));
data_interp=datetime(TSG_INTERP_data, 'ConvertFrom', 'datenum');

figure
set(gcf,'color','w');
scatter(data_interp(d),ANOM_jpl(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
c=colorbar;
cmocean('balance')
%caxis([-5 6])
% colormap(redbluecmap)
c.Label.String = ('Anomalies (ºC)');
caxis([-5 5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Anomalies (ºC)')
title('Temperature Anomalies')
subtitle('Southern Ocean - 2015', 'FontSize', 14);
grid on
%xline(data, '--r')

figure
set(gcf,'color','w');
scatter(data_interp(d),TSG_INTERP(d), 25, ANOM_jpl(d), 'filled', 'MarkerEdgeColor','k')
ylim([4 22])
c=colorbar;
cmocean('balance')
c.Label.String = ('Anomalies (ºC)');
caxis([-5 5])
set(gca,'fontsize',18)
xlabel('Time')
ylabel('Temperature (ºC)')
title('Temperature Anomalies')
subtitle('Southern Ocean - 2015', 'FontSize', 14);
grid on
%xline