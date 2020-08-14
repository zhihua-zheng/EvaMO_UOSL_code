function phi_lstar_E6(pzData,iD,E6,axm,fpart)

%% Setting

zet_range(1,:) = [2e-3 4];
zet_range(2,:) = [3e-2 4];
phi_range(1,:) = [6e-2 1];
phi_range(2,:) = [6e-2 1];

lstar_range = [0.9 1.5];

bin_lim(1,:) = [-3.0 0.7];
bin_lim(2,:) = [-1.6 0.7];

mkrNum{1} = '(a)';
mkrNum{2} = '(b)';
mkrNum{3} = '(c)';
mkrNum{4} = '(d)';
mkrAlign{1} = 'left';
mkrAlign{2} = 'right';
mkrLoc = [0.02 0.98];

hcb_pos = [0.105 0.095 0.33 0.025];

off_r = 0.02;
iL = iD + 2;
cAlpha = 0.35;
cStr{1} = 'orangish';
cStr{2} = 'easter purple';
cStr{3} = 'aqua blue';

%% Unpack data

phi  = pzData.phi_qs;
zeta = pzData.zeta_qs;
etaX = pzData.etaX_qs;
etaY = pzData.etaY_qs;
fzS  = pzData.fzS_qs;

%% Preprocessing

Igd = ~isnan(phi) & ~isnan(zeta) & ~isnan(fzS) & ~isnan(etaX) &...
      ~isnan(etaY);

x_data = zeta(Igd);
y_data = phi(Igd);
f_data = fzS(Igd);
etX_data = etaX(Igd);
etY_data = etaY(Igd);

Inx = x_data <= 0;

nx_data = -x_data(Inx);
ny_data =  y_data(Inx);
nf_data =  f_data(Inx);
nex_data = etX_data(Inx);
ney_data = etY_data(Inx);

%% Computation

% "forced convection"
Ifrcc = nx_data <= 3;

% Empirical dimensionless functions
Nzet        = -logspace(-5,1.5,100);
[phiS71,~] = get_emp_phi(Nzet,'Kansas');

[phiS,~,Nlstar] = get_H15se_phi(-nx_data,nex_data,ney_data,nf_data,E6,2);
[phiSo,~] = get_H15se_phi(-nx_data,nex_data,ney_data,nf_data,6,2);

% average length scale ratio l/k|z| in bins
Betxi = -1:0.1:0.2;
Bzeti = -2:0.2:0;
Betxc = (Betxi(1:end-1) + Betxi(2:end))/2;
Bzetc = (Bzeti(1:end-1) + Bzeti(2:end))/2;

nZet = length(Bzetc);
nEtx = length(Betxc);
Blc  = zeros(nEtx,nZet);

[~,~,~,bZet,bEtx] = histcounts2(-nx_data,nex_data,Bzeti,Betxi);

for j = 1:nZet
    for i = 1:nEtx
        ind = and(bZet == j, bEtx == i);
        if sum(~isnan(Nlstar(ind))) >= 5
            Blc(i,j) = nanmean(Nlstar(ind));
        end
    end
end

%% Figure

switch fpart  
    case 1 % Unstable side

% -------------------------------------------------------------------------
% ----- phi plot --------------------------------------------------------
% -------------------------------------------------------------------------
bin_xi = logspace(bin_lim(iD,1),bin_lim(iD,2),20)';
Ngood  = numel(nx_data(Ifrcc));
NminC  = Ngood*0.005;

Nref71 = plot(axm(iD),-Nzet,phiS71,'color',[.5 .5 .5],...
                      'linewidth',3);
hold(axm(iD),'on')
grid(axm(iD),'on')

% Langmuir turbulence model with E6 = 6
Szet = pinBin(bin_xi,nx_data(Ifrcc),nx_data(Ifrcc),'left');
Sphi = pinBin(bin_xi,nx_data(Ifrcc),phiSo(Ifrcc),'left');
Sphi.qm(Sphi.n < NminC) = NaN;
nLTo = plot(axm(iD),Szet.qm,Sphi.qm,'d','MarkerSize',15,'LineWidth',1,...
      'MarkerFaceColor',rgb(cStr{1}),'MarkerEdgeColor',rgb('orange'));
errorbar(axm(iD),Szet.qm,Sphi.qm,Sphi.qmL,Sphi.qmU,'d','CapSize',5,...
      'MarkerEdgeColor',rgb('orange'),'MarkerSize',15,...
      'Color',rgb('orange'),'linewidth',2.5,'linestyle','none');

% Langmuir turbulence model with E6 = 2.5
Sphi = pinBin(bin_xi,nx_data(Ifrcc),phiS(Ifrcc),'left');
Sphi.qm(Sphi.n < NminC) = NaN;
nLT = plot(axm(iD),Szet.qm,Sphi.qm,'d','MarkerSize',15,'LineWidth',1,...
      'MarkerFaceColor',rgb(cStr{2}),'MarkerEdgeColor',rgb('dark lavender'));
errorbar(axm(iD),Szet.qm,Sphi.qm,Sphi.qmL,Sphi.qmU,'d','CapSize',5,...
      'MarkerEdgeColor',rgb('dark lavender'),'MarkerSize',15,...
      'Color',rgb('dark lavender'),'linewidth',2.5,'linestyle','none');

% observations
Sphi = pinBin(bin_xi,nx_data(Ifrcc),ny_data(Ifrcc),'left');
Sphi.qm(Sphi.n < NminC) = NaN;
nOBS = plot(axm(iD),Szet.qm,Sphi.qm,'d','MarkerSize',15,'LineWidth',1,...
       'MarkerFaceColor',rgb(cStr{3}),'MarkerEdgeColor',rgb('sea blue'));
errorbar(axm(iD),Szet.qm,Sphi.qm,Sphi.qmL,Sphi.qmU,'d','CapSize',5,...
      'MarkerEdgeColor',rgb('sea blue'),'MarkerSize',15,...
      'Color',rgb('sea blue'),'linewidth',2.5,'linestyle','none');

text(axm(iD),mkrLoc(iD),0.98,mkrNum{iD},'Units','Normalized','FontSize',22,...
     'HorizontalAlignment',mkrAlign{iD},'VerticalAlignment','top')
set(axm(iD),'Xscale','log','Yscale','log','Xdir','reverse',...
    'TickDir','out','Box','on',...
    'xlim',zet_range(iD,:),'ylim',phi_range(iD,:))

% adjust label position
axm(iD).Position = axm(iD).Position + off_r*[1 1 -1 -1];
xlh = xlabel(axm(iD),['-',char(950),' = -|z|/L'],'fontsize',18);
ylh = ylabel(axm(iD),[char(981),'_h'],'fontsize',18,'Interpreter','tex');
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
set(xlh,'Position',xlh.Position + [0 -off_r 0]);
set(ylh,'Position',ylh.Position + [-off_r 0 0]);

% transparent markers
setMarkerColor(nLTo,cStr{1},cAlpha);
setMarkerColor(nLT, cStr{2},cAlpha);
setMarkerColor(nOBS,cStr{3},cAlpha);

% legend
if iD == 2
    nlgd = legend(axm(iD),[Nref71 nLTo nLT nOBS],...
   'Prediction: Monin-Obukhov (Kansas curve)',...
   'SE Langmuir turbulence model & length scale Eq. (E_6=6)',...
   'SE Langmuir turbulence model & length scale Eq. (E_6=2.5)',...
   'Observation: gradient from fitting',...
   'FontSize',16,'AutoUpdate','off','Interpreter','tex',...
   'Position',[0.25 0.835 0.61 0.13]);

    set(nlgd.BoxFace,'ColorType','truecoloralpha',...
        'ColorData',uint8(255*[1; 1; 1; .8]))
    legendMarkers([nLTo nLT nOBS],nlgd,cStr,cAlpha)
    
    set(axm(iD),'YTickLabels',[])
    delete(ylh)
end

% -------------------------------------------------------------------------
% ----- lstar plot --------------------------------------------------------
% -------------------------------------------------------------------------
cmap_seq = cbrewer('seq','GnBu',6);
histogram2(axm(iL),'XBinEdges',Bzeti,'YBinEdges',Betxi,...
               'BinCounts',Blc','DisplayStyle','tile','FaceAlpha',.8,...
               'EdgeColor','k','EdgeAlpha',.5,'LineWidth',1);
colormap(axm(iL),cmap_seq); caxis(axm(iL),lstar_range)
if iD == 1
    hcb = colorbar(axm(iL),'northoutside','AxisLocation','in',...
                   'FontSize',14,'LineWidth',1.5);
    hcb.Position = hcb_pos;
    clh = ylabel(hcb,[char(8467),'/',char(954),'|z|'],'fontsize',18);
    clh.Units = 'Normalized';
    clh.HorizontalAlignment = 'center';
    clh.Position = [0.5 0 0];
end

hold(axm(iL),'on')
plot(axm(iL),[-2 0],[0 0],'--k')
text(axm(iL),mkrLoc(iD),0.98,mkrNum{iL},'Units','Normalized','FontSize',22,...
     'HorizontalAlignment',mkrAlign{iD},'VerticalAlignment','top')
text(axm(iL),0.02,0.66,'E_6=2.5','Units','Normalized','FontSize',16,...
     'HorizontalAlignment','left','VerticalAlignment','bottom',...
     'Interpreter','tex','Color',rgb('liliac'),'FontWeight','bold')
set(axm(iL),'TickDir','out','Box','on','xlim',[-2 0],'ylim',[-1 0.2])

% adjust label position
axm(iL).Position = axm(iL).Position + off_r*[1 1 -1 -1];
xlh = xlabel(axm(iL),[char(950),' = |z|/L'],'fontsize',18); 
ylh = ylabel(axm(iL),['Normalized downwind Stokes shear ',char(951),...
                      '^x'],'fontsize',18,'Interpreter','tex');
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
set(xlh,'Position',xlh.Position + [0 -off_r 0]);
set(ylh,'Position',ylh.Position + [-off_r 0 0]);

if iD == 2
    set(axm(iL),'YTickLabels',[])
    delete(ylh)
end

end