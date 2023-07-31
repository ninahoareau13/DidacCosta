function [] = plot_SOFRESH_drake(x,y,z,tt)

figure 
%set(gcf,'position',[249 408 780 567],'color','w');
set(gcf,'position',[249 608 1400 1200],'color','w');
%m_proj('miller','lat',[-77 77]);
%m_proj('azimuthal equal-area','latitude',-60,'longitude',-60,'radius',20,'rectbox','on');

%m_proj('lambert','long',[9 31],'lat',[53 67]);
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(x,y,z);shading flat;hold on;colormap(jet);
%m_grid('xtick',16,'tickdir','out','ytick',[50 60 70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
c = colorbar('fontsize',30,'location','southoutside');
c.Label.String =(tt);
set(gca,'fontsize',30)
m_grid('box','fancy','linestyle','-','gridcolor','k');

%m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[.2 .65 1]);

%m_grid('xticklabel',[],'linewi',2,'box','fancy','tickdir','in','Fontsize',24);
hold on