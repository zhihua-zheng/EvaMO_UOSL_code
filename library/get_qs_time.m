%% get_qs_time

%% Compute time scales

% free convection velocity scale
SKF.Wstar = nthroot(-BLD .* BfH',3);

Ub = max([SKF.Wstar SKF.Ustar],[],2,'includenan');
Tturnover = BLD./Ub/3600; % [hour]

itime_mid = (itime(1:end-1) + itime(2:end))/2;

% buoyancy forcing time scale
dBo_dt = center_diff(BfH',itime,1,'mid');
dBo_dtl = [nan; dBo_dt];
dBo_dtr = [dBo_dt; nan];
TBol = abs(BfH' ./ dBo_dtl)*24; % [hour]
TBor = abs(BfH' ./ dBo_dtr)*24; % [hour]
TBo  = mean([TBol TBor],2,'includenan');

% time scale of boundary layer evolution
dD_dt = center_diff(BLD,itime,1,'mid'); % [m/hr]
dD_dtl = [nan; dD_dt];
dD_dtr = [dD_dt; nan];
Tdl = abs(BLD ./ dD_dtl)*24; % [hour]
Tdr = abs(BLD ./ dD_dtr)*24; % [hour]

% time scale of wind stress variation
dwstr_dt = center_diff(SKF.Ustar.^2,itime,1,'mid');
dwstr_dtl = [nan; dwstr_dt];
dwstr_dtr = [dwstr_dt; nan];
Twstrl = abs(SKF.Ustar.^2 ./ dwstr_dtl)*24; % [hour]
Twstrr = abs(SKF.Ustar.^2 ./ dwstr_dtr)*24; % [hour]
Twstr  = mean([Twstrl Twstrr],2,'includenan');

% external forcing time scale
Texf = min([TBo Twstr],[],2,'includenan');

%% Quasi-steady stages

stage = nan(size(itime));
stage(Texf > 10*Tturnover) = 1;

%% Show quasi-steady state periods

if iD == 2

figure('position',[0 0 900 600]);
[hax,~] = tight_subplot(2,1,[.03 .02],[.1 .05],[.1 .1]);

j = 3;

nj = j*250;
pj = (j-1)*250 + 1;
Jchunk = pj:nj;
Jqs = Jchunk(stage(Jchunk) == 1);

% upper panel
axes(hax(1))
yyaxis left
plot(idatm_loc(Jchunk),-BfH(Jchunk),'-k','linewidth',1.5); hold on
plot(idatm_loc([pj,nj]),[0 0],'--k')
ylim([-9e-7 3e-7]); xlim(idatm_loc([pj,nj])); grid on
text(0.02,0.04,'(a)','Units','Normalized','FontSize',22,...
     'HorizontalAlignment','left','VerticalAlignment','bottom')
set(hax(1),'XTickLabels',[])
hax(1).YAxis(1).Color = 'k';
ylh1 = ylabel('-B_f^H [m^2 s^{-3}]','fontsize',18,'Interpreter','tex');
ylh1.Units = 'Inches';
set(ylh1,'Position',ylh1.Position - [0.2 0 0]);

yyaxis right
plot(idatm_loc(Jchunk),1025*SKF.Ustar(Jchunk).^2,'linewidth',1.5,...
     'color',rgb('purple')); hold on
ylim([0 0.3])
for jj = 1:length(Jqs)
    plot_sq_patch(idatm_loc(Jqs(jj)),[0 0.3])
end
hax(1).YAxis(2).Color = rgb('purple');
ylh2 = ylabel([char(964),'_w',' [N m^{-2}]'],...
              'fontsize',18,'Interpreter','tex');
ylh2.Units = 'Inches';
set(ylh2,'Position',ylh2.Position + [0.2 0 0]);


% lower panel
axes(hax(2))
plot(idatm_loc(Jchunk),Texf(Jchunk),':k','linewidth',2); hold on
plot(idatm_loc(Jchunk),10*Tturnover(Jchunk),...
     'color',rgb('azure'),'linewidth',1.5)
for jj = 1:length(Jqs)
    plot_sq_patch(idatm_loc(Jqs(jj)),[0 80])
end
xtickformat('dd'); grid on
text(0.02,0.96,'(b)','Units','Normalized','FontSize',22,...
     'HorizontalAlignment','left','VerticalAlignment','top')
ylim([0 80]); xlim(idatm_loc([pj,nj]))

ylh3 = ylabel('time scale [hr]','fontsize',18);
ylh3.Units = 'Inches';
set(ylh3,'Position',ylh3.Position - [0.3 0 0]);
legend('external forcing','boundary layer eddy * 10','fontsize',16,...
       'location','north')
   
set(hax,'TickDir','out')
end