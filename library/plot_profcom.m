function plot_profcom(rpData,iD,axm)

%% Setting

off_r = 0.02;
cAlpha = 0.2;

mkrNum{1} = '(b)';
mkrNum{2} = '(d)';

%% Unpack data

rObs  = rpData.rObs;  
rMOi  = rpData.rMOi;
zbl   = rpData.zbl;
Iprof = rpData.Iprof;

%% Preprocessing

tmp1 = zbl(:,Iprof);
tmp2 = rObs(:,Iprof);
tmp3 = rMOi(:,Iprof);

bin_xi = [(-1:.04:-0.32) (-0.29:.03:-0.14) (-0.12:.02:-0.04) ...
          (-0.03:.015:0)]';
Sobs = pinBin(bin_xi,tmp1(:),tmp2(:),'left');
Szml = pinBin(bin_xi,tmp1(:),tmp1(:),'left');
Smoi = pinBin(bin_xi,tmp1(:),tmp3(:),'left');

Sobs.qm(Sobs.n < 3) = NaN;
Smoi.qm(Smoi.n < 3) = NaN;

% referenced to mean temperature at SLD
Jref = find(Szml.qm<-0.2,1,'last');
tmp2 = tmp2 - Sobs.qm(Jref);
Sobs.qm = Sobs.qm-Sobs.qm(Jref);
tmp3 = tmp3 - Smoi.qm(Jref);
Smoi.qm = Smoi.qm-Smoi.qm(Jref);

f_ox = [Sobs.qm-Sobs.qmL; flip(Sobs.qm+Sobs.qmU)];
f_mx = [Smoi.qm-Smoi.qmL; flip(Smoi.qm+Smoi.qmU)];
f_y  = [Szml.qm; flip(Szml.qm)];

%% Figure

fill(axm,1e-1*[-1 -1 1 1],[0 -0.2 -0.2 0],rgb('squash'),...
     'LineStyle','none','FaceAlpha',.05); 
hold(axm,'on')
grid(axm,'on')
plot(axm,[0 0],[-1 0],'--','color',.7*ones(1,3),'linewidth',2)     

% plot all profiles
scatter(axm,tmp2(:),tmp1(:),2,'MarkerFaceColor','none',...
        'MarkerEdgeColor',rgb('lightblue'),'MarkerEdgeAlpha',cAlpha);
scatter(axm,tmp3(:),tmp1(:),2,'MarkerFaceColor','none',...
        'MarkerEdgeColor',rgb('soft pink'),'MarkerEdgeAlpha',cAlpha);

% plot ensemble averages
plot(axm,Sobs.qm(~isnan(Sobs.qm)),Szml.qm(~isnan(Sobs.qm)),'color',...
         rgb('azure'),'linewidth',3)
fill(axm,f_ox(~isnan(f_ox)),f_y(~isnan(f_ox)),...
         rgb('azure'),'LineStyle','none','FaceAlpha',.3)
plot(axm,Smoi.qm(~isnan(Smoi.qm)),Szml.qm(~isnan(Smoi.qm)),'color',...
         rgb('coral'),'linewidth',3)
fill(axm,f_mx(~isnan(f_mx)),f_y(~isnan(f_mx)),...
         rgb('coral'),'LineStyle','none','FaceAlpha',.3)

set(axm,'ylim',[-0.27 0],'xlim',[-2.2e-2 .5e-2],'TickDir','out',...
        'YAxisLocation','right')

text(axm,0.98,0.98,mkrNum{iD},'Units','Normalized','FontSize',22,...
     'HorizontalAlignment','right','VerticalAlignment','top')

ylh = ylabel(axm,'z/H','fontsize',18);
xlh = xlabel(axm,['Temp - Temp_{SL} [',char(0176),'C]'],'fontsize',18);
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
axm.Position = axm.Position + off_r*[0 1 -1 -1];
set(xlh,'Position',xlh.Position + [0 -off_r 0]);
set(ylh,'Position',ylh.Position + [off_r 0 0]); 

if iD == 1    
    set(axm,'XTickLabels',[])
    delete(xlh)
end

end