function PZ_scat_obs(phi_obs,zet_obs)

%% Load data

zet_range = [3e-3 8];       
phi_range = [2e-2 2.0];

x1 = zet_obs{1};
x2 = zet_obs{2};
y1 = phi_obs{1};
y2 = phi_obs{2};
y1(y1<0) = NaN;
y2(y2<0) = NaN;

Igd1 = ~isnan(x1) & ~isnan(y1);
Igd2 = ~isnan(x2) & ~isnan(y2);

x1_data = x1(Igd1);
y1_data = y1(Igd1);
x2_data = x2(Igd2);
y2_data = y2(Igd2);

Inx1 = x1_data <= 0;
Inx2 = x2_data <= 0;

nx1_data = -x1_data(Inx1);
ny1_data =  y1_data(Inx1);
nx2_data = -x2_data(Inx2);
ny2_data =  y2_data(Inx2);

%% Computation

% Empirical dimensionless functions
Nzet        = -logspace(-5,1.5,100);
[Nphi_71,~] = get_emp_phi(Nzet,'Kansas');

nx_edge = (-2.4:0.4:0.8);
Num_nx  = numel(nx_edge);
nxc     = (nx_edge(1:end-1) + nx_edge(2:end))/2;

ny1_grp = nan(size(ny1_data));
[~,~,bin1] = histcounts(log10(nx1_data),nx_edge);
for i = 1:length(nxc); ny1_grp(bin1==i) = nxc(i); end

ny2_grp = nan(size(ny2_data));
[~,~,bin2] = histcounts(log10(nx2_data),nx_edge);
for i = 1:length(nxc); ny2_grp(bin2==i) = nxc(i); end

Ibox1 = ismember(nxc,ny1_grp);
Ibox2 = ismember(nxc,ny2_grp);

%% Figure

cStr{1} = 'azure';
cStr{2} = 'coral';
off_ratio = 0.02;

figure('position',[0 0 570 500]);
Naxm = axes;
Nref_71 = plot(Naxm,log10(-Nzet),log10(Nphi_71),...
               'color',[.5 .5 .5],'linewidth',3); hold on

boxplot(Naxm,log10(ny1_data),ny1_grp,'Positions',nxc(Ibox1)+0.05,...
        'Widths',0.09,'Notch','on');
set(Naxm,'XTickLabel',{' '},'XTick',[])
boxplot(Naxm,log10(ny2_data),ny2_grp,'Positions',nxc(Ibox2)-0.05,...
        'Widths',0.09,'Notch','on');
set(Naxm,'XTickLabel',{' '},'XTick',[])

boxG = findobj(Naxm,'Tag','boxplot'); % order is reversed
rebox_prop(boxG(2),log10(ny1_data),ny1_grp,nxc(Ibox1),0.99,cStr{1});
rebox_prop(boxG(1),log10(ny2_data),ny2_grp,nxc(Ibox2),0.99,cStr{2});

y_tick = log10([1e-2 1e-1 1]);
y_ticl = {'10^{-2}','10^{-1}','10^{0}'};
yticks(y_tick); yticklabels(y_ticl); ylim(log10(phi_range))
Naxm.YAxis.MinorTick = 'on';
Naxm.YAxis.MinorTickValues = log10([1e-2*(2:9) 1e-1*(2:9) (2:9)]);
x_tick = log10([1e-3 1e-2 1e-1 1 1e1]);
x_ticl = {'10^{-3}','10^{-2}','10^{-1}','10^{0}'};
xticks(x_tick); xticklabels(x_ticl); xlim(log10(zet_range))
Naxm.XAxis.MinorTick = 'on';
Naxm.XAxis.MinorTickValues = log10([1e-3*(2:9) 1e-2*(2:9) 1e-1*(2:9)...
                                    (2:9) 1e1*(2:9)]);
plot(Naxm,log10([1e-5  5e1]),[0 0],'--k')
plot(Naxm,[nx_edge; nx_edge],log10([2e-2 2.3e-2]'*ones(1,Num_nx)),....
    'linewidth',2,'color','k')
grid(Naxm,'on')
set(Naxm,'xdir','reverse','TickDir','out','TickLabelInterpreter','tex',...
         'xminorgrid','on','yminorgrid','on','Box','on')

boxph = findobj(Naxm,'Type','Patch');
% sort by xposition
[~,ind] = sort(cellfun(@mean, get(boxph, 'XData')),'descend');
nlgd = legend(Naxm,[Nref_71 boxph(ind(1:2))'],...
         'Prediction: Monin-Obukhov (Kansas curve)',...
         'Observation: OCSP',...
         'Observation: SPURS-I',...
         'fontsize',16,'location','northwest','AutoUpdate','off');

set(nlgd.BoxFace, 'ColorType','truecoloralpha',...
    'ColorData',uint8(255*[1; 1; 1; .8]))

% adjust label positions
xlh = xlabel(Naxm,['-',char(950),' = -|z|/L'],'fontsize',18);
ylh = ylabel(Naxm,[char(981),'_h'],'fontsize',18,'Interpreter','tex');
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
Naxm.Position = Naxm.Position + off_ratio*[1 1 -1 -1];
set(xlh,'Position',xlh.Position + [0 -off_ratio 0]);
set(ylh,'Position',ylh.Position + [-off_ratio 0 0]);

end