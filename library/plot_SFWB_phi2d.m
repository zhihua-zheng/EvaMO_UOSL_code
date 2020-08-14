function plot_SFWB_phi2d(axm,r2o,pzPA,pzSP)

%% Load observational data

zetaPA = pzPA.zeta_qs;
xiPA   = pzPA.xi_qs;
xiPA(xiPA<1) = NaN;
Icvc   = ~isnan(pzPA.phi_qs) & ~isnan(pzPA.fzS_qs) & ...
         ~isnan(pzPA.etaX_qs) & ~isnan(xiPA) & zetaPA<0;
zetaPA = zetaPA(Icvc);
xiPA   = xiPA(Icvc);

zetaSP = pzSP.zeta_qs;
xiSP   = pzSP.xi_qs;
xiSP(xiSP<1) = NaN;
Icvc   = ~isnan(pzSP.phi_qs) & ~isnan(pzSP.fzS_qs) & ...
         ~isnan(pzSP.etaX_qs) & ~isnan(xiSP) & zetaSP<0;
zetaSP = zetaSP(Icvc);
xiSP   = xiSP(Icvc);

%% Computation

zeta = -2:0.02:0;
xi   =  1:.1:10;
[Zeta,Xi] = meshgrid(zeta,xi);

zeta_range = [min(zeta) max(zeta)];
xi_range   = [min(xi)   max(xi)];

[pS0,~] = get_SMCse_phi_nol(zeta,'KC1994');
[pS1,~] = get_SFWB_phi_apr(Zeta,Xi,1);
% [pS2,~]   = get_SFWB_phi_de(Zeta,Xi);

%% Figure for Phi_h, with approximation

pS_lev = [.1 .22 .25 .3 .45 .6]';
cmap_seq = cbrewer('seq','Purples',10);

[~,h1] = contourf(axm,zeta,xi,pS1,pS_lev,'linewidth',2,...
                  'linecolor',[.75 .75 .75]);
hold(axm,'on')
colormap(axm,cmap_seq); caxis(axm,[0 1])
hcb = colorbar(axm,'northoutside');
hcb.Position = [0.073 0.91 0.46 0.03];

hsPA = scatter(axm,zetaPA,xiPA,10,rgb('azure'),'s','filled',...
               'MarkerEdgeColor',[.6 .6 .6],'MarkerFaceAlpha',0.25,...
               'MarkerEdgeAlpha',0.2);
hsSP = scatter(axm,zetaSP,xiSP,10,rgb('coral'),'s','filled',...
               'MarkerEdgeColor',[.6 .6 .6],'MarkerFaceAlpha',0.25,...
               'MarkerEdgeAlpha',0.2);
plot(axm,r2o*[1 10],xi_range,'color',rgb('amber'),'linewidth',3)

contour(axm,zeta,xi,pS1,pS_lev,'linewidth',2,'linecolor',[.75 .75 .75])
[c2,h2] = contour(axm,zeta,xi,pS1./pS0,[.2 .4 .6 .8 1 1.2 1.4],...
                  'linecolor','k','showtext','on','linewidth',1.2);
clabel(c2,h2,'fontsize',14,'fontname','SansSerif','LabelSpacing',505)

grid(axm,'on')
set(axm,'TickDir','out','Box','on')
text(axm,0.02,0.02,'(a)','Units','Normalized','FontSize',22,...
         'HorizontalAlignment','left','VerticalAlignment','bottom')
set(axm,'ydir','reverse','xlim',zeta_range,'ylim',xi_range)
hl1 = legend(axm,[h1 h2 hsPA hsSP],...
       [char(981),'_h(',char(950),', ',char(958),')'],...
       [char(981),'_h(',char(950),', ',char(958),')/',char(981),...
        '_h^{no wave}'],...
       'OCSP data','SPURS-I data',...
       'fontsize',16,'Position',[0.08 0.2 0.19 0.2],'Interpreter','tex');
set(hl1.BoxFace,'ColorType','truecoloralpha', ...
    'ColorData',uint8(255*[1; 1; 1; .6]))

% modify scatter marker in legend
drawnow;
hl1_comp = hl1.EntryContainer.Children;
set(hl1_comp(1).Icon.Transform.Children.Children,...
    'Size',8, 'FaceColorData',uint8(255*[rgb('coral')'; 0.7]))
set(hl1_comp(2).Icon.Transform.Children.Children,...
    'Size',8, 'FaceColorData',uint8(255*[rgb('azure')'; 0.7]))

% adjust label position
off_r = 0.02;
axm.Position = axm.Position + off_r*[1 1.1 -0.9 -2];
xlh = xlabel(axm,[char(950),' = |z|/L'],  'fontsize',18,'Interpreter','tex');
ylh = ylabel(axm,[char(958),' = |z|/z_0'],'fontsize',18,'Interpreter','tex');
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
set(xlh,'Position',xlh.Position + [0 -off_r 0]);

%% Figure for Phi_h, no approximation

% pS_lev = [0 .22 .25 .3 .45 .6]';
% cmap_seq = cbrewer('seq','Purples',10);
% 
% figure('position',[0 0 530 480]);
% [~,h1] = contourf(zeta,xi,pS2,pS_lev,'linewidth',1.2,...
%                   'linecolor',[.6 .6 .6]); hold on
% colormap(cmap_seq); caxis([0 1]); colorbar
% % colorbar('Ticks',pS_lev,...
% %          'TickLabels',cellstr(num2str(round(pS_lev,2))))
% 
% hsPA = scatter(zetaPA,xiPA,10,rgb('azure'),'s','filled',...
%                'MarkerFaceAlpha',0.2);
% hsSP = scatter(zetaSP,xiSP,10,rgb('coral'),'s','filled',...
%                'MarkerFaceAlpha',0.2);
% 
% contour(zeta,xi,pS2,pS_lev,'linewidth',1.2,'linecolor',[.6 .6 .6])
% [c2,h2] = contour(zeta,xi,pS2./pS0,[.2 .4 .6 .8 1 1.1 1.2 1.3],...
%                   'linecolor','k','showtext','on','linewidth',1.2);
% clabel(c2,h2,'fontsize',13,'fontname','SansSerif','LabelSpacing',490)
% 
% % plot([0 0],xi_range,':k')
% box on; grid on
% set(gca,'gridlinestyle',':','gridalpha',.3,'fontsize',12,'TickDir','out')
% text(0.02,0.98,'(a)','Units','Normalized','FontSize',22,...
%      'HorizontalAlignment','left','VerticalAlignment','top')
% xlim(zeta_range); ylim(xi_range)
% hl2 = legend([h1 h2 hsPA hsSP],...
%        [char(981),'_h(',char(950),', ',char(958),')'],...
%        [char(981),'_h(',char(950),', ',char(958),')/',char(981),...
%         '_h^{no wave}'],...
%        'OCSP data','SPURS-I data',...
%        'fontsize',16,'Position',[0.15 0.6 0.32 0.24],'Interpreter','tex');
% set(hl2.BoxFace, 'ColorType','truecoloralpha', ...
%     'ColorData',uint8(255*[1; 1; 1; .8])) % 0.8 means 20% transparent
% 
% drawnow;
% hl2_comp = hl2.EntryContainer.Children;
% set(hl2_comp(1).Icon.Transform.Children.Children,...
%     'Size',8, 'FaceColorData',uint8(255*[rgb('coral')'; 0.7]))
% set(hl2_comp(2).Icon.Transform.Children.Children,...
%     'Size',8, 'FaceColorData',uint8(255*[rgb('azure')'; 0.7]))
% 
% % adjust label position
% axm = gca;
% axm.Position = axm.Position + off_r*[1 1 -1 -1];
% xlh = xlabel(char(950),'fontsize',18,'Interpreter','tex');
% ylh = ylabel(char(958),'fontsize',18,'Interpreter','tex');
% xlh.Units = 'Normalized';
% ylh.Units = 'Normalized';
% set(xlh,'Position',xlh.Position + [0 -off_r 0]);
% set(ylh,'Position',ylh.Position + [-off_r 0 0]);

%% Figure for Phi_m

% pM_lev = [0 .2 .4 .6 .8 1 1.2 1.4]';
% cmap_seq = flip( cbrewer('div','RdYlGn',10) );
% 
% figure('position',[0 0 530 480]);
% [~,h1] = contourf(zeta,xi,pM1,pM_lev,'linewidth',1,...
%                   'linecolor',[.9 .9 .9]); hold on
% colormap(cmap_seq); caxis([0 2]); colorbar
% % colorbar('Ticks',pM_lev,...
% %          'TickLabels',cellstr(num2str(pM_lev)))
% [c2,h2] = contour(zeta,xi,pM1./pM0,[.2 .4 .6 .8 1],...
%                   'linecolor','k','showtext','on');
% clabel(c2,h2,'fontsize',13,'fontname','SansSerif','LabelSpacing',450)
% % plot([0 0],-xi_range,':k')
% box on; grid on
% set(gca,'gridlinestyle',':','gridalpha',.3,'TickDir','out')
% xlim(zeta_range); ylim(xi_range)
% hl3 = legend([h1 h2],...
%        [char(981),'_m(',char(950),', ',char(958),')'],...
%        [char(981),'_hm',char(950),', ',char(958),')/',char(981),...
%         '_h^{no wave}'],...
%        'fontsize',16,'Position',[0.15 0.6 0.34 0.26],'Interpreter','tex');
% set(hl3.BoxFace, 'ColorType','truecoloralpha', ...
%     'ColorData',uint8(255*[1; 1; 1; .8])) % 0.8 means 20% transparent
% 
% % adjust label position
% axm = gca;
% axm.Position = axm.Position + off_r*[1 1 -1 -1];
% xlh = xlabel(char(950),'fontsize',18,'Interpreter','tex');
% ylh = ylabel(char(958),'fontsize',18,'Interpreter','tex');
% xlh.Units = 'Normalized';
% ylh.Units = 'Normalized';
% set(xlh,'Position',xlh.Position + [0 -off_r 0]);
% set(ylh,'Position',ylh.Position + [-off_r 0 0]);

end