function demo_sfwbProf(axm,r2o)

%% Constants

BLD = 60;
SLD = BLD/5;
Hs  = 2;
wspd = 10;
kappa = 0.4;
cp = 3985;
alpha = 1.66e-4;

%% Surface wave breaking model

z0   = Hs*0.6;
zdum = -SLD:.1:-z0;
xi   = abs(zdum)/z0;

MOL = z0/r2o;
Ustar = sqrt(stresstc(wspd,10)/1025);
Bo   = Ustar^3/kappa/MOL;
zeta = abs(zdum)/MOL;

Qo = Bo/9.81/alpha;
heatF = Qo*1025*cp;
Tstar = Qo/Ustar;

pSsfwb = get_SFWB_phi_apr(zeta,xi,1);
pSMO   = get_emp_phi(zeta,'Kansas');

intgrd = Tstar/kappa ./ abs(zdum) .* pSsfwb;
sfwbT = cumtrapz(zdum,intgrd);

intgrd = Tstar/kappa ./ abs(zdum) .* pSMO;
MOT = cumtrapz(zdum,intgrd);

%% Figure: MOST predicted prof vs SFWB model predicted prof

p1 = plot(axm,MOT,xi,'color',[.5 .5 .5],'linewidth',4);
p1.Color(4) = 0.7;
grid(axm,'on'); hold(axm,'on')
p2 = plot(axm,sfwbT,xi,'linewidth',4,'color',rgb('amber'));
p2.Color(4) = 0.7;

set(axm,'ydir','reverse','Box','on','YAxisLocation','right',...
        'XAxisLocation','Top','Color','none','YTickLabel',[],...
        'Xlim',[-0.08 0],'TickDir','out')

end