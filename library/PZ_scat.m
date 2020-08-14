function PZ_scat(pzData,iD,fh,axm,fpart)

%% Setting

zet_range(1,:) = [2e-3 5];
zet_range(2,:) = [3e-2 5];
zet_range(3,:) = [3e-3 5];

bin_lim(1,:) = [-3.0 0.7];
bin_lim(2,:) = [-1.6 0.7];
bin_lim(3,:) = [-3.0 0.7];

phi_range(1,:) = [6e-2 1.5];
phi_range(2,:) = [6e-2 1.5];
phi_range(3,:) = [6e-2 1.5];

nBin = 4;
mkrNum{1} = '(a)';
mkrNum{2} = '(b)';
mkrAlign{1} = 'left';
mkrAlign{2} = 'right';
mkrLoc = [0.02 0.98];

off_r = 0.02;
cAlpha = 0.35;
cStr{1} = 'amber';
cStr{2} = 'amber';
cStr{3} = 'salmon';
cStr{4} = 'aqua blue';

%% Unpack data

phi  = pzData.phi_qs;  
zeta = pzData.zeta_qs;
etaX = pzData.etaX_qs;
etaY = pzData.etaY_qs;
fzS  = pzData.fzS_qs;
xi   = pzData.xi_qs;

%% Preprocessing

xi(xi<1) = NaN;
Igood = ~isnan(phi) & ~isnan(zeta) & ~isnan(etaX) &...
        ~isnan(etaY) & ~isnan(fzS) & ~isnan(xi);

x_data   = zeta(Igood);
y_data   = phi(Igood);
etX_data = etaX(Igood);
etY_data = etaY(Igood);
fzS_data = fzS(Igood);
xi_data  = xi(Igood);

IposX = x_data >  0;
InegX = x_data <= 0;

py_data  = y_data(IposX);
px_data  = x_data(IposX);
pex_data = etX_data(IposX);
pey_data = etY_data(IposX);
pf_data  = fzS_data(IposX);
pxi_data = xi_data(IposX);

ny_data = y_data(InegX);
nx_data = -x_data(InegX);
nex_data = etX_data(InegX);
ney_data = etY_data(InegX);
nf_data = fzS_data(InegX);
nxi_data = xi_data(InegX);

%% Figure

switch fpart
    case 1 % Unstable side

Ifrcc  = nx_data <= 3; % "forced convection"
bin_xi = logspace(bin_lim(iD,1),bin_lim(iD,2),24)';
Ngood  = numel(nx_data(Ifrcc));
NminC  = Ngood*0.005;

% Empirical dimensionless functions
Nzet       = -logspace(-5,1.5,100);
[Nphi71,~] = get_emp_phi(Nzet,'Kansas');

Nphi_H15  = get_H15se_phi_nol(-nx_data,nex_data,ney_data,nf_data,2);
Nphi_H15B = get_H15se_phi_nol(-nx_data,nex_data,ney_data,0,2);
Nphi_H15T = get_H15se_phi_nol(-nx_data,nex_data,ney_data,1,2);
Nphi_WBv1 = get_SFWB_phi_apr(-nx_data,nxi_data,1);

Nref71 = plot(axm,-Nzet,Nphi71,'color',[.5 .5 .5],...
                  'linewidth',3);
hold(axm,'on')
grid(axm,'on')
scatter(axm,nx_data,ny_data,1,'filled',...
        'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',.1,...
        'MarkerEdgeColor',[.5 .5 .5],'MarkerEdgeAlpha',.3);

% Langmuir turbulence model predictions with constant fzs
SphiH15B = pinBin(bin_xi,nx_data(Ifrcc),Nphi_H15B(Ifrcc),'left');
SphiH15T = pinBin(bin_xi,nx_data(Ifrcc),Nphi_H15T(Ifrcc),'left');
Szet     = pinBin(bin_xi,nx_data(Ifrcc),nx_data(Ifrcc),'left');
SphiH15B.qm(SphiH15B.n < NminC) = NaN;
SphiH15T.qm(SphiH15T.n < NminC) = NaN;

Ib = ~isnan(SphiH15B.qmL) & ~isnan(SphiH15B.qmU) & ~isnan(SphiH15B.qm);
It = ~isnan(SphiH15T.qmL) & ~isnan(SphiH15T.qmU) & ~isnan(SphiH15T.qm);

plot(axm,Szet.qm,SphiH15B.qm,'color',rgb('rosa'),'linewidth',1)
fill(axm,[Szet.qm(Ib);flip(Szet.qm(Ib))],...
    [SphiH15B.qm(Ib)-SphiH15B.qmL(Ib);flip(SphiH15B.qm(Ib)+SphiH15B.qmU(Ib))],...
     rgb('rosa'),'LineStyle','none','FaceAlpha',cAlpha)
plot(axm,Szet.qm,SphiH15T.qm,'color',rgb('rosa'),'linewidth',1)
fill(axm,[Szet.qm(It);flip(Szet.qm(It))],...
    [SphiH15T.qm(It)-SphiH15T.qmL(It);flip(SphiH15T.qm(It)+SphiH15T.qmU(It))],...
     rgb('rosa'),'LineStyle','none','FaceAlpha',cAlpha)

% wave breaking model prediction
Sphi = pinBin(bin_xi,nx_data(Ifrcc),Nphi_WBv1(Ifrcc),'left');
Szet = pinBin(bin_xi,nx_data(Ifrcc),nx_data(Ifrcc),  'left');
Sphi.qm(Sphi.n < NminC) = NaN;
nWB1 = plot(axm,Szet.qm,Sphi.qm,'d','MarkerSize',15,'LineWidth',1,...
       'MarkerFaceColor',rgb(cStr{1}),'MarkerEdgeColor',rgb('dirty orange'));
errorbar(axm,Szet.qm,Sphi.qm,Sphi.qmL,Sphi.qmU,'d','CapSize',5,...
         'MarkerEdgeColor',rgb('dirty orange'),'MarkerSize',15,...
         'Color',rgb('dirty orange'),'linewidth',2.5,'linestyle','none');

% Langmuir turbulence model prediction
Sphi = pinBin(bin_xi,nx_data(Ifrcc),Nphi_H15(Ifrcc),'left');
Sphi.qm(Sphi.n < NminC) = NaN;
nLT = plot(axm,Szet.qm,Sphi.qm,'d','MarkerSize',15,'LineWidth',1,...
      'MarkerFaceColor',rgb(cStr{3}),'MarkerEdgeColor',rgb('faded red'));
errorbar(axm,Szet.qm,Sphi.qm,Sphi.qmL,Sphi.qmU,'d','CapSize',5,...
         'MarkerEdgeColor',rgb('faded red'),'MarkerSize',15,...
         'Color',rgb('faded red'),'linewidth',2.5,'linestyle','none');

% observations
Sphi = pinBin(bin_xi,nx_data(Ifrcc),ny_data(Ifrcc),'left');
Sphi.qm(Sphi.n < NminC) = NaN;
nOBS = plot(axm,Szet.qm,Sphi.qm,'d','MarkerSize',15,'LineWidth',1,...
       'MarkerFaceColor',rgb(cStr{4}),'MarkerEdgeColor',rgb('sea blue'));
errorbar(axm,Szet.qm,Sphi.qm,Sphi.qmL,Sphi.qmU,'d','CapSize',5,...
         'MarkerEdgeColor',rgb('sea blue'),'MarkerSize',15,...
         'Color',rgb('sea blue'),'linewidth',2.5,'linestyle','none');

plot(axm,[1e-5  5e1],[1 1],'--k')
text(axm,mkrLoc(iD),0.98,mkrNum{iD},'Units','Normalized','FontSize',22,...
     'HorizontalAlignment',mkrAlign{iD},'VerticalAlignment','top')
set(axm,'Xscale','log','Yscale','log','Xdir','reverse',...
        'TickDir','out','Box','on')

% adjust label position
axm.Position = axm.Position + off_r*[1 1 -1 -1];
xlh = xlabel(axm,['-',char(950),' = -|z|/L'],'fontsize',18);
ylh = ylabel(axm,[char(981),'_h'],'fontsize',18,'Interpreter','tex');
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
set(xlh,'Position',xlh.Position + [0 -off_r 0]);
set(ylh,'Position',ylh.Position + [-off_r 0 0]);

% PDF on x-axis
axmP    = axm.Position;
ax_xpdf = axes('Parent',fh,'Position',[axmP(1:3) axmP(4)/7]);
plot_xpdf_log(ax_xpdf,nx_data,'gray')
linkaxes([axm ax_xpdf],'x')

% PDF on y-axis
ax_ypdf = axes('Parent',fh,'Position',[axmP(1:2) axmP(3)/7 axmP(4)]);
plot_ypdf_log(ax_ypdf,ny_data,'gray')
linkaxes([axm ax_ypdf],'y')

% xlim and ylim after adding PDF
xlim(axm,zet_range(iD,:)); ylim(axm,phi_range(iD,:))

% transparent markers
setMarkerColor(nWB1,cStr{1},cAlpha);
setMarkerColor(nLT, cStr{3},cAlpha);
setMarkerColor(nOBS,cStr{4},cAlpha);

% legend
if iD == 2
    nlgd = legend(axm,[Nref71 nWB1 nLT nOBS],...
   'Prediction: Monin-Obukhov (Kansas curve)',...
  ['"SE" surface wave breaking model (',char(8467),'=',char(954),'|z|)'],...
  [' SE  Langmuir turbulence model (',char(8467),'=',char(954),'|z|)'],...
   'Observation: gradient from fitting',...
   'FontSize',16,'AutoUpdate','off',...
   'Position',[0.4 0.73 0.38 0.2]);

    set(nlgd.BoxFace,'ColorType','truecoloralpha',...
        'ColorData',uint8(255*[1; 1; 1; .8]))
    legendMarkers([nWB1 nLT nOBS],nlgd,cStr([1,3:4]),cAlpha)
    
    set(axm,'YTickLabels',[])
    delete(ylh)
end


    case 2 % Stable side
    
Pzet       = logspace(-5,1.5,100);
[Pphi71,~] = get_emp_phi(Pzet,'Kansas');
Pphi_H15  = get_H15se_phi_nol(px_data,pex_data,pey_data,pf_data,2);
Pphi_H15B = get_H15se_phi_nol(px_data,pex_data,pey_data,0,2);
Pphi_H15T = get_H15se_phi_nol(px_data,pex_data,pey_data,1,2);

Pgood = numel(px_data);
PminC = max(Pgood*0.001,5);

bin_xi = logspace(bin_lim(iD,1),bin_lim(iD,2),30)';
Szet   = pinBin(bin_xi,px_data,px_data,'left');

Pref71 = plot(axm,Pzet,Pphi71,'-k', 'linewidth',3); 
hold(axm,'on')
grid(axm,'on')
scatter(axm,px_data,py_data,2,'filled',...
        'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',.1,...
        'MarkerEdgeColor',[.5 .5 .5],'MarkerEdgeAlpha',.3);        

SphiH15B = pinBin(bin_xi,px_data,Pphi_H15B,'left');
SphiH15T = pinBin(bin_xi,px_data,Pphi_H15T,'left');
SphiH15B.qm(SphiH15B.n < PminC) = NaN;
SphiH15T.qm(SphiH15T.n < PminC) = NaN;

plot(axm,Szet.qm,SphiH15B.qm,'color',rgb('rosa'),'linewidth',1.5)
plot(axm,Szet.qm,SphiH15T.qm,'color',rgb('rosa'),'linewidth',1.5)

% Langmuir turbulence model
Sphi = pinBin(bin_xi,px_data,Pphi_H15,'left');
Sphi.qm(Sphi.n < PminC) = NaN;
pLT = plot(axm,Szet.qm,Sphi.qm,'d','MarkerSize',15,'LineWidth',1.5,...
      'MarkerFaceColor',rgb(cStr{3}),'MarkerEdgeColor',rgb('faded red'));
setMarkerColor(pLT,cStr{3},cAlpha);
errorbar(axm,Szet.qm,Sphi.qm,Sphi.qmL,Sphi.qmU,'d','CapSize',10,...
         'MarkerEdgeColor',rgb('faded red'),'MarkerSize',15,...
         'Color',rgb('faded red'),'linewidth',1.5,'linestyle','none');

% observational data
Sphi = pinBin(bin_xi,px_data,py_data,'left');
if sum(Sphi.n >= PminC) >= nBin
    Sphi.qm(Sphi.n < PminC) = NaN;
    
    pOBS = plot(axm,Szet.qm,Sphi.qm,'d','MarkerSize',15,'LineWidth',1.5,...
           'MarkerFaceColor',rgb(cStr{4}),'MarkerEdgeColor',rgb('sea blue'));
    setMarkerColor(pOBS,cStr{4},cAlpha);
    errorbar(axm,Szet.qm,Sphi.qm,Sphi.qmL,Sphi.qmU,'d','CapSize',10,...
             'MarkerEdgeColor',rgb('sea blue'),'MarkerSize',15,...
             'Color',rgb('sea blue'),'linewidth',1.5,'linestyle','none');
else
    close(fp)
    return
end

plot(axm,[1e-5  5e1],[1 1],'--k')
text(axm,mkrLoc(iD),0.98,mkrNum{iD},'Units','Normalized','FontSize',22,...
     'HorizontalAlignment',mkrAlign{iD},'VerticalAlignment','top')
set(axm,'Xscale','log','Yscale','log','TickDir','out','Box','on')

% adjust label position
axm.Position = axm.Position + off_r*[1 1 -1 -1];
xlh = xlabel(axm,[char(950),' = |z|/L'],'fontsize',18);
ylh = ylabel(axm,[char(981),'_h'],'fontsize',18,'Interpreter','tex');
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
set(xlh,'Position',xlh.Position + [0 -off_r 0]);
set(ylh,'Position',ylh.Position + [-off_r 0 0]);

xlim(axm,zet_range(iD,:)); ylim(axm,[6e-2 1e1])

if iD == 2
    plgd = legend(axm,[Pref71 pLT pOBS],...
           'Prediction: Monin-Obukhov (Kansas curve)',...
           'Prediction: SMC with Langmuir turbulence',...
           'Observation: gradient from fitting',...
           'FontSize',16,...
           'Position',[0.3 0.8 0.45 0.15]);

    set(plgd.BoxFace,'ColorType','truecoloralpha',...
        'ColorData',uint8(255*[1; 1; 1; .8]))
    legendMarkers([pLT pOBS],plgd,cStr(3:4),cAlpha)
    
    set(axm,'YTickLabels',[])
    delete(ylh)
end

end
%% Subfunction for plotting horizontal pdf on log scale

function plot_xpdf_log(axPdf,pdfData,cStr)

histf = figure;
Shis  = histogram(log10(pdfData),'Normalization','pdf');
pdf   = Shis.Values;
npdf  = pdf/max(pdf);
bi    = Shis.BinEdges;
b     = (bi(1:end-1)+bi(2:end))/2;

plot(axPdf,10.^b,npdf,'color',rgb(cStr))
hold(axPdf,'on')
fill(axPdf,[10.^b flip(10.^b)],[zeros(1,Shis.NumBins) flip(npdf)],...
     rgb(cStr),'LineStyle','none','FaceAlpha',.15)
set(axPdf,'Xscale','log','xdir','reverse',...
    'XTick',[],'YTick',[],'Color','none','Visible','off')

close(histf)
end

%% Subfunction for plotting vertical pdf on log scale

function plot_ypdf_log(axPdf,pdfData,cStr)

pdfData(pdfData<0) = NaN;
histf = figure;
Shis  = histogram(log10(pdfData),'Normalization','pdf');
pdf   = Shis.Values;
npdf  = pdf/max(pdf);
bi    = Shis.BinEdges;
b     = (bi(1:end-1)+bi(2:end))/2;

plot(axPdf,npdf,10.^b,'color',rgb(cStr))
hold(axPdf,'on')
fill(axPdf,[zeros(1,Shis.NumBins) flip(npdf)],[10.^b flip(10.^b)],...
     rgb(cStr),'LineStyle','none','FaceAlpha',.15)
set(axPdf,'Yscale','log','XTick',[],'YTick',[],...
    'Color','none','Visible','off')

close(histf)
end

end