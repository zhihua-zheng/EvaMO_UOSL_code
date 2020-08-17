function plot_H15se_phi2d(LTax,mop,pzPA,pzSP)

%% Load observational data

zetaPA = pzPA.zeta_qs;
etaxPA = pzPA.etaX_qs;
fzsPA  = pzPA.fzS_qs;
Icvc   = ~isnan(pzPA.phi_qs) & ~isnan(fzsPA) & ~isnan(pzPA.xi_qs) &...
         ~isnan(etaxPA) & zetaPA<0;
zetaPA = zetaPA(Icvc);
etaxPA = etaxPA(Icvc);
fzsPA  = fzsPA(Icvc);

zetaSP = pzSP.zeta_qs;
etaxSP = pzSP.etaX_qs;
fzsSP  = pzSP.fzS_qs;
Icvc   = ~isnan(pzSP.phi_qs) & ~isnan(fzsSP) & ~isnan(pzSP.xi_qs) &...
         ~isnan(etaxSP) & zetaSP<0;
zetaSP = zetaSP(Icvc);
etaxSP = etaxSP(Icvc);
fzsSP  = fzsSP(Icvc);

%% Computation

zet = -2:0.02:0;
etX = -1:0.01:0.2;
fzs = (0:0.02:1);

cEtY = 0;
cEtX = -0.5;
cfzS = 0.3;

[zeta,etaX] = meshgrid(zet,etX);
zet_range = [min(zet) max(zet)];
eta_range = [min(etX) max(etX)];

if strcmp(mop,'H2013')
    [pS0,~] = get_H13se_phi_nol(zet,0,0,2);
    [pS1,~] = get_H13se_phi_nol(zeta,etaX,cEtY,2);

elseif strcmp(mop,'H2015')
    [pS0,~] = get_H15se_phi_nol(zet,0,0,0,2);
    [pS1,~] = get_H15se_phi_nol(zeta,etaX,cEtY,cfzS,2);

else
    disp('function type not available!');
    return
end

%% Figure: Phi_h in (zeta, eta^x) space

pS_lev = [.1 .15 .2 .3 .4 .5 .7]';
cmap_seq = cbrewer('seq','Purples',10);

[~,h1] = contourf(LTax(1),zet,etX,pS1,pS_lev,'linewidth',2,...
                  'linecolor',[.7 .7 .7]);
hold(LTax(1),'on')
colormap(LTax(1),cmap_seq); caxis(LTax(1),[0 1]); 
hcb1 = colorbar(LTax(1),'northoutside');
hcb1.Position = [0.082 0.92 0.39 0.03];

hsPA = scatter(LTax(1),zetaPA,etaxPA,10,rgb('azure'),'s','filled',...
               'MarkerEdgeColor',[.6 .6 .6],'MarkerFaceAlpha',0.25,...
               'MarkerEdgeAlpha',0.2);          
hsSP = scatter(LTax(1),zetaSP,etaxSP,10,rgb('coral'),'s','filled',...
               'MarkerEdgeColor',[.6 .6 .6],'MarkerFaceAlpha',0.25,...
               'MarkerEdgeAlpha',0.2);          

contour(LTax(1),zet,etX,pS1,pS_lev,'linewidth',2,'linecolor',[.7 .7 .7]);
[c2,h2] = contour(LTax(1),zet,etX,pS1./pS0,[.5 .6 .7 .8 .9 1 1.1 1.2],...
                  'linecolor','k','showtext','on','linewidth',1.2);
clabel(c2,h2,'fontsize',14,'fontname','SansSerif','LabelSpacing',370)

plot(LTax(1),zet_range,[0 0],':k')
plot(LTax(1),zet_range,cEtX*[1 1],'--','color',rgb('amber'),'linewidth',3)
grid(LTax(1),'on')
set(LTax(1),'TickDir','out','Box','on','xlim',zet_range,'ylim',eta_range)
text(LTax(1),0.02,0.98,'(a)','Units','Normalized','FontSize',22,...
         'HorizontalAlignment','left','VerticalAlignment','top')

hl1 = legend([h1 h2 hsPA hsSP],...
      [char(981),'_h(',char(950),', ',char(951),'^x) with ',...
       char(951),'^y = 0, ',char(402),'^{s}_z = ',num2str(cfzS)],...
      [char(981),'_h(',char(950),', ',char(951),'^x)/',...
       char(981),'_h^{no wave}'],...
      'OCSP data','SPURS-I data',...
      'fontsize',16,'location','southwest','Interpreter','tex');
set(hl1.BoxFace,'ColorType','truecoloralpha', ...
    'ColorData',uint8(255*[1; 1; 1; .7]))
    
% modify scatter marker in legend
drawnow;
hl1_comp = hl1.EntryContainer.Children;
set(hl1_comp(1).Icon.Transform.Children.Children,...
    'Size',8, 'FaceColorData',uint8(255*[rgb('coral')'; 0.7]))
set(hl1_comp(2).Icon.Transform.Children.Children,...
    'Size',8, 'FaceColorData',uint8(255*[rgb('azure')'; 0.7]))

% adjust label position
off_r = 0.02;
LTax(1).Position = LTax(1).Position + off_r*[1 1 -1 -1];
xlh = xlabel(LTax(1),[char(950),' = |z|/L'],'fontsize',18);
ylh = ylabel(LTax(1),['Normalized downwind Stokes shear ',char(951),...
                      '^x'],'fontsize',18,'Interpreter','tex');
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
set(xlh,'Position',xlh.Position + [0 -off_r 0]);
set(ylh,'Position',ylh.Position + [-off_r 0 0]);

%% Figure for Phi_m

% pM_lev = [-.2 -.1 0 .1 .2 .3 .4 .5 .75 1 1.5]';
% cmap_seq = flip( cbrewer('div','RdYlGn',11) );
% 
% figure('position',[0 0 520 480]);
% [~,h1] = contourf(zet,etX,pM1,pM_lev,'linewidth',1,...
%                   'linecolor',[.9 .9 .9]); hold on
% colormap(cmap_seq); caxis([-.2 2]); colorbar
% % colorbar('Ticks',pM_lev,...
% %          'TickLabels',cellstr(num2str(pM_lev)))
% [c2,h2] = contour(zet,etX,pM1./pM0,[-.2 0 .2 .4 .6 .8 1 2],...
%                   'linecolor','k','showtext','on');
% clabel(c2,h2,'fontsize',13,'fontname','Courier','LabelSpacing',450)
% plot(zet_range,[0 0],':k');plot([0 0],eta_range,':k')
% box on; grid on
% set(gca,'gridlinestyle',':','gridalpha',.3,'fontsize',12,'TickDir','out')
% xlim(zet_range); ylim(eta_range)
% xlabel('$\zeta$', 'fontsize',18,'Interpreter','latex')
% ylabel('$\eta^x$','fontsize',18,'Interpreter','latex')
% hl2 = legend([h1 h2],'$\phi_m (\zeta, \eta^x)$',...
%                '$\phi_m (\zeta, \eta^x) / \phi_m (\zeta)$',...
%                'fontsize',20,'location','best','Interpreter','latex');
% set(hl2.BoxFace, 'ColorType','truecoloralpha', ...
%     'ColorData',uint8(255*[1; 1; 1; .75])) % 0.8 means 20% transparent
% title(['$\phi_m$ in ',textModel],...
%        'fontsize',20,'Interpreter','latex')

%% Figure for length scale ratio rs

% rs_lev = [0.2 0.4 0.6 .8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8]';
% cmap_seq = cbrewer('seq','GnBu',15);
% 
% figure('position',[0 0 520 480]);
% [c1,h1] = contourf(zet,etX,rs1,rs_lev,'linewidth',1.2,...
%                   'linecolor',[.4 .4 .4],'showtext','on'); hold on
% colormap(cmap_seq); 
% caxis([0 3]); colorbar
% clabel(c1,h1,'fontsize',19,'fontname','Courier','fontangle','italic',...
%       'LabelSpacing',450)
% 
% plot(zet_range,[0 0],':k')
% box on; grid on
% set(gca,'TickDir','out')
% xlim(zet_range); ylim(eta_range)
% xlabel('$\zeta$', 'fontsize',18,'Interpreter','latex')
% ylabel('$\eta^x$','fontsize',18,'Interpreter','latex')
% hl3 = legend(h1,'$l / \kappa |z|$',...
%        'fontsize',20,'location','best','Interpreter','latex');
% set(hl3.BoxFace, 'ColorType','truecoloralpha', ...
%     'ColorData',uint8(255*[1; 1; 1; .75])) % 0.8 means 20% transparent

% title(['length scale ratio in ',textModel],...
%        'fontsize',20,'Interpreter','latex')

%% Variation of Phi_h in (zeta, fzs) space

[zeta,fzS] = meshgrid(zet,fzs);
[pS2,~] = get_H15se_phi_nol(zeta,cEtX,0,fzS,2);

%% Figure: phi_h in (zeta, fzs) space

pS_lev = [.1 .15 .2 .3 .4 .5 .7]';
cmap_seq = cbrewer('seq','Purples',10);

[~,h1] = contourf(zet,fzs,pS2,pS_lev,'linewidth',2,...
                  'linecolor',[.7 .7 .7]); 
hold(LTax(2),'on')
colormap(LTax(2),cmap_seq); caxis(LTax(2),[0 1]);
hcb2 = colorbar(LTax(2),'northoutside');
hcb2.Position = [0.507 0.92 0.39 0.03];

hsPA = scatter(zetaPA,fzsPA,10,rgb('azure'),'s','filled',...
               'MarkerEdgeColor',[.6 .6 .6],'MarkerFaceAlpha',0.25,...
               'MarkerEdgeAlpha',0.2);          
hsSP = scatter(zetaSP,fzsSP,10,rgb('coral'),'s','filled',...
               'MarkerEdgeColor',[.6 .6 .6],'MarkerFaceAlpha',0.25,...
               'MarkerEdgeAlpha',0.2);          

contour(zet,fzs,pS2,pS_lev,'linewidth',2,'linecolor',[.7 .7 .7])
[c2,h2] = contour(zet,fzs,pS2./pS0,[.5 .6 .7 .8 .9 1 1.1 1.2],...
                  'linecolor','k','showtext','on','linewidth',1.2);
clabel(c2,h2,'fontsize',14,'fontname','SansSerif','LabelSpacing',435)

plot(LTax(2),zet_range,cfzS*[1 1],'--','color',rgb('amber'),'linewidth',3)
grid(LTax(2),'on')
set(LTax(2),'TickDir','out','Box','on','xlim',zet_range,'ylim',[0 1],...
            'YAxisLocation','right')
text(LTax(2),0.02,0.98,'(b)','Units','Normalized','FontSize',22,...
         'HorizontalAlignment','left','VerticalAlignment','top')

hl4 = legend([h1 h2 hsPA hsSP],...
      [char(981),'_h(',char(950),', ',char(402),'^{s}_z) with ',...
       char(951),'^y = 0, ',char(951),'^x = ',num2str(cEtX)],...
      [char(981),'_h(',char(950),', ',char(402),'^{s}_z)/',...
       char(981),'_h^{no wave}'],...
      'OCSP data','SPURS-I data',...
      'fontsize',16,'location','northwest','Interpreter','tex');
hl4.Position = [hl4.Position(1) hl4.Position(2)-0.11 hl4.Position(3:4)];
set(hl4.BoxFace,'ColorType','truecoloralpha', ...
    'ColorData',uint8(255*[1; 1; 1; .7]))

% modify scatter marker in legend
drawnow;
hl4_comp = hl4.EntryContainer.Children;
set(hl4_comp(1).Icon.Transform.Children.Children,...
    'Size',8, 'FaceColorData',uint8(255*[rgb('coral')'; 0.7]))
set(hl4_comp(2).Icon.Transform.Children.Children,...
    'Size',8, 'FaceColorData',uint8(255*[rgb('azure')'; 0.7]))

% adjust label position
off_r = 0.02;
LTax(2).Position = LTax(2).Position + off_r*[0 1 -1 -1];
xlh = xlabel(LTax(2),[char(950),' = |z|/L'],'fontsize',18);
ylh = ylabel(LTax(2),['Surface proximity function ',char(402),...
                      '^{s}_z'],'fontsize',18,'Interpreter','tex');
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
set(xlh,'Position',xlh.Position + [0 -off_r 0]);
set(ylh,'Position',ylh.Position + [off_r 0 0]);

end