%% demo_spursi_cprof

%% Plot profiles

j = Ict(319);
f_prof = figure('position',[ 0 0 480 690]);
ax1 = axes('Parent',f_prof);

% salinity
Ig  = ~isnan(iSAprof(:,j));
sp1 = plot(ax1,iSAprof(Ig,j),-depth_s(Ig),'s-',...
               'color',rgb('green teal'),'linewidth',2);
sp1.Color(4) = 0.3;
hold on
sp2 = plot(ax1,cSA_ict(Ig,319),-depth_s(Ig),'-.',...
               'color',rgb('gray'),'linewidth',3);
sp2.Color(4) = 0.7;
sp3 = plot(ax1,SAprof(Ig,j),-depth_s(Ig),...
               'color',rgb('green teal'),'linewidth',4);
sp3.Color(4) = 0.5;

% temperature
Ig = ~isnan(iPTprof(:,j));
tp1 = plot(ax1,1.45*iPTprof(Ig,j),-depth_t(Ig),'s-',...
               'color',rgb('azure'),'linewidth',2);
tp1.Color(4) = 0.3;
tp3 = plot(ax1,1.45*PTprof(Ig,j),-depth_t(Ig),...
               'color',rgb('azure'),'linewidth',4);
tp3.Color(4) = 0.5;

% density
Ig = ~isnan(iPDprof(:,j));
rp1 = plot(ax1,1.505*iPDprof(Ig,j),-depth_t(Ig),'s-',...
               'color',rgb('amethyst'),'linewidth',2);
rp1.Color(4) = 0.3;
rp3 = plot(ax1,1.505*PDprof(Ig,j),-depth_t(Ig),...
               'color',rgb('amethyst'),'linewidth',4);
rp3.Color(4) = 0.5;

%% Adjust axis properties

ax1.GridColorMode = 'manual';
grid(ax1,'on')
ax1.GridAlpha = 0.3;
ax1.XColor = rgb('green teal');
axPos = ax1.Position;
ax1.Position = axPos + [0 0.22 0 -0.22];

% cover the top XAxis with an empty axis
ax0 = axes('Parent',f_prof,'Position',ax1.Position,... 
      'XAxisLocation','top','YAxisLocation','left',...
      'XTickLabels',[],'YTickLabels',[],'Color','none');
ax0.XTick = ax1.XTick;
ax0.YTick = ax1.YTick;
ax0.XAxis.TickValuesMode = 'auto';
ax0.YAxis.TickValuesMode = 'auto';
linkaxes([ax0 ax1])

%% Add axes for temperatuer and density

ax2 = axes('Parent',f_prof,'position',(axPos.*[1 1 1 1e-3])+[0 .11 0 0],...
           'color','none','XColor',rgb('azure'),'YColor',rgb('azure'));
ax3 = axes('Parent',f_prof,'position',(axPos.*[1 1 1 1e-3])+[0  0  0 0],...
           'color','none','XColor',rgb('amethyst'),'YColor',rgb('amethyst'));
set([ax0 ax1 ax2 ax3],'fontsize',14,'TickDir','out','Units','Normalized')

ax1.YLim = [-60 0];
ax1.XLim = [37.35 37.95];
ax2.XLim = [37.35 37.95]/1.45;
ax3.XLim = [37.35 37.95]/1.505;
xticks(ax3,[24.85 24.95 25.05 25.15])

ybl1 = ylabel(ax1,'z [m]','FontSize',18);
xbl1 = xlabel(ax1,'Salinity [g kg^{-1}]','FontSize',18,'Interpreter','tex');
xbl2 = xlabel(ax2,['Temperature [',char(0176),'C]'],'FontSize',18);
xbl3 = xlabel(ax3,'Density [kg m^{-3}]','FontSize',18,'Interpreter','tex');

ybl1.Units = 'Normalized';
xbl1.Units = 'Inches';
xbl2.Units = 'Inches';
xbl3.Units = 'Inches';

off_r = 0.02;
ax0.Position = ax0.Position + off_r*[1 0 -1  0];
ax1.Position = ax1.Position + off_r*[1 0 -1  0];
ax2.Position = ax2.Position + off_r*[1 0 -1  0];
ax3.Position = ax3.Position + off_r*[1 0 -1  0];

set(ybl1,'Position',ybl1.Position + off_r*[-1 0 0])
set(xbl1,'Position',xbl1.Position - [0 .1 0])
set(xbl2,'Position',xbl2.Position - [0 .1 0])
set(xbl3,'Position',xbl3.Position - [0 .1 0])
