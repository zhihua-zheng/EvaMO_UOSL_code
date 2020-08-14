function slp = dT_scat(dTData,Islc,iD,axm,fpart,Fitc)

%% Setting

modeStr = 'Observations';
mkrNum{1} = '(a)';
mkrNum{2} = '(c)';
cR(1,:) = [0 0.70 0.93 1];
cR(2,:) = [0 0.20 0.93 1];
fb(1,:) = [-0.022 0.005];
fb(2,:) = [-0.022 0.005];

tbound = fb(iD,:);
xltext = ['MOST predictions [',char(0176),'C]'];
yltext = [modeStr,' [',char(0176),'C]'];
mycolormap = customcolormap(cR(iD,:),...
             [rgb('bluish');    rgb('azure');...
              rgb('pale blue'); rgb('heather')]);
off_r = 0.02;

% option to include intercept in linear fit
if Fitc
    itc_opt = true;
else
    itc_opt = false;
end

%% Unpack data

dT_MOi = dTData.dT_MOi;
dT_fit = dTData.dT_fit;
dT_obs = dTData.dT_obs;

%% Preprocessing

dT_MOi_qs = dT_MOi(:,Islc);
dT_fit_qs = dT_fit(:,Islc);

Iunst    = dT_MOi    <= 0;
Istab    = dT_MOi    >  0;
Iunst_qs = dT_MOi_qs <= 0;
Istab_qs = dT_MOi_qs >  0;

unstX   = dT_MOi(Iunst);
unstY   = dT_obs(Iunst);
unstXqs = dT_MOi_qs(Iunst_qs);
unstYqs = dT_fit_qs(Iunst_qs);

stabX   = dT_MOi(Istab);
stabY   = dT_obs(Istab);
stabXqs = dT_MOi_qs(Istab_qs);
stabYqs = dT_fit_qs(Istab_qs);

unstXqs = unstXqs(:);
unstYqs = unstYqs(:);
unstX   = unstX(:);
unstY   = unstY(:);

stabXqs = stabXqs(:);
stabYqs = stabYqs(:);
stabX   = stabX(:);
stabY   = stabY(:);

%% Figure

switch fpart
    case 1 % Unstable side

Xedges = (-3.3e-2:1.5e-4:6e-3);
Yedges = (-3.3e-2:1.5e-4:6e-3);
[binVal,~,~] = histcounts2(unstXqs,unstYqs,Xedges,Yedges,...
                           'Normalization','count');
binVal = binVal/max(binVal(:));
histogram2(axm,'XBinEdges',Xedges,'YBinEdges',Yedges,...
           'BinCounts',binVal,'DisplayStyle','tile','FaceAlpha',.6,...
           'EdgeColor','none');
       
hold(axm,'on')

% percentage of dT_obs distribution for each contour line
vl  = [0.3 0.6 0.9 0.96];
nvl = numel(vl);

[Hist,Xi,Yi,~,~] = histcounts2(unstX,unstY,'Normalization','count');
Xctr = (Xi(1:end-1) + Xi(2:end))/2;
Yctr = (Yi(1:end-1) + Yi(2:end))/2;
Hsum = sum(Hist(:));
Hper = sort(Hist(:),'descend') ./ Hsum;
Hcum = cumsum(Hper);

% find the count value corrsponding to the percentage of distribution in vl
hvl = nan(1,nvl);
for i = 1:nvl
    [~,ind] = min(abs(Hcum-vl(i)));
    hvl(i)  = Hper(ind);
end

prob = Hist./Hsum;
prob(prob==0) = 1e-12;
contour(axm,Xctr,Yctr,prob',hvl,'linecolor',0.5*ones(1,3),'linewidth',1.6)

np(1) = sum(~isnan(unstYqs));
cbin_min = 0.005*np(1);
nBin     = 3;

if np(1) >= cbin_min*nBin
    % bin average
    [cbin,bin_xi] = histcounts(unstXqs(unstXqs >= tbound(1)),...
                               'BinWidth',0.0015);

    if sum(cbin >= cbin_min) < nBin
        [cbin,bin_xi] = histcounts(unstXqs(unstXqs >= tbound(1)),...
                                   'BinWidth',0.001);
        if sum(cbin >= cbin_min) < nBin
            [~,bin_xi] = histcounts(unstXqs(unstXqs >= tbound(1)),...
                                    'BinWidth',0.0005);
        end
    end

    nSy = pinBin(bin_xi',unstXqs,unstYqs,'left');
    nSx = pinBin(bin_xi',unstXqs,unstXqs,'left');

    if sum(nSy.n >= cbin_min & ~isnan(nSy.qm)) >= nBin
        nSy.qm(nSy.n < cbin_min) = NaN;
        
        % unstable_side_fit
        xref = linspace(tbound(1),0);
        lfm  = fitlm(nSx.qm,nSy.qm,'linear','RobustOpts','on',...
                     'Intercept',itc_opt);
        if Fitc
            incpt = lfm.Coefficients.Estimate(1);
            slp   = lfm.Coefficients.Estimate(2);
        else
            incpt = 0;
            slp   = lfm.Coefficients.Estimate(1);
        end
        [yref,yCI] = predict(lfm,xref');
        
        nfitStr = ['slope = ',num2str(round(slp,2),'% 4.2f'),...
                 ', intercept = ',num2str(round(incpt,4),'% 6.4f')];
        fill(axm,[xref flip(xref)]',[yCI(:,1);flip(yCI(:,2))],...
                 rgb('coral'),'LineStyle','none','FaceAlpha',.3)
        nfit = plot(axm,xref,yref,'Color',rgb('coral'),'lineWidth',1.8);
        scatter(axm,nSx.qm,nSy.qm,200,'filled',...
                    'MarkerFaceColor',rgb('yellow'))
        errorbar(axm,nSx.qm,nSy.qm,nSy.qmL,nSy.qmU,'o','MarkerSize',13,...
                 'CapSize',10,'Color',rgb('ultramarine'),...
                 'linewidth',1.2,'linestyle','none')
    end
end

xlim(axm,tbound); ylim(axm,tbound)
plot(axm,[-.1 .2],[-.1 .2],'--','color',.7*ones(1,3),'linewidth',2)
plot(axm,[-.1 .2],[0 0],':k'); plot(axm,[0 0],[-.1 .2],':k')
set(axm,'TickDir','out','Box','on')
nlgd = legend(axm,nfit,nfitStr,'location','southwest',...
              'fontsize',16,'AutoUpdate','off');
set(nlgd.BoxFace,'ColorType','truecoloralpha',...
    'ColorData',uint8(255*[1; 1; 1; .8]))

% adjust label position
axm.YAxis.Exponent = -3;
axm.XAxis.Exponent = -3;
xlh = xlabel(axm,xltext,'fontsize',18);
ylh = ylabel(axm,yltext,'fontsize',18);
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
axm.Position = axm.Position + off_r*[1 1 -1 -1];
set(xlh,'Position',xlh.Position + [0 -off_r 0]);
set(ylh,'Position',ylh.Position + [-off_r 0 0]);

text(axm,0.02,0.98,mkrNum{iD},'Units','Normalized','FontSize',22,...
     'HorizontalAlignment','left','VerticalAlignment','top')

% colorbar location
colormap(axm,mycolormap)
hc1 = colorbar(axm,'LineWidth',1.5,...
                   'Position',[axm.Position(1)+axm.Position(3)-0.075 ...
                               axm.Position(2)+0.047 0.03 0.33]);
clh = ylabel(hc1,'Normalized count','fontsize',18,'color','w');
clh.Units = 'Normalized';
clh.HorizontalAlignment = 'center';
clh.Position = [0 0.6 0];
caxis(axm,[0 0.1])

if iD == 1    
    set(axm,'XTickLabels',[])
    delete(xlh)
end

    % PDF for MOST prediction
%     Nax_pos = Naxm.Position;
%     Nax_xpdf = axes('Position',[Nax_pos(1:3) Nax_pos(4)/2]);
%     plot_xpdf(Nax_xpdf,unstX(~isnan(unstX)),'gray')
%     plot_xpdf(Nax_xpdf,unstY(~isnan(unstY)),'azure')
%     linkaxes([Naxm,Nax_xpdf],'x')

    % PDF for obs
%     Nax_ypdf = axes('Position',[Nax_pos(1:2) Nax_pos(3)/7 Nax_pos(4)]);
%     plot_ypdf(Nax_ypdf,unstY,'azure')
%     linkaxes([Naxm,Nax_ypdf],'y')

    case 2 % Stabel side

Xedges = [-Inf -5e-3:1e-4:1e-1 Inf];
Yedges = [-Inf -5e-3:1e-4:1e-1 Inf];
[binVal,~,~] = histcounts2(stabXqs,stabYqs,Xedges,Yedges,...
                           'Normalization','count');
binVal = binVal/max(binVal(:));
    
np(2) = sum(~isnan(stabYqs));
cbin_min = 0.005*np(2);
nBin     = 3;

if np(2) >= max(cbin_min*nBin,30)
    
    histogram2(axm,'XBinEdges',Xedges,'YBinEdges',Yedges,...
               'BinCounts',binVal,'DisplayStyle','tile','FaceAlpha',.6,...
               'EdgeColor','none'); hold on
    colormap(mycolormap); hc2 = colorbar; caxis([0 1])
    ylabel(hc2,'Normalized count','fontsize',18)

    % bin average
    [cbin,bin_xi] = histcounts(stabXqs(stabXqs <= 1e-1),...
                               'BinWidth',0.0015);

    if sum(cbin >= cbin_min) < nBin
        [cbin,bin_xi] = histcounts(stabXqs(stabXqs <= 1e-1),...
                                   'BinWidth',0.001);
        if sum(cbin >= cbin_min) < nBin
            [~,bin_xi] = histcounts(stabXqs(stabXqs <= 1e-1),...
                                    'BinWidth',0.0005);
        end
    end

    pSy = pinBin(bin_xi',stabXqs,stabYqs,'left');
    pSx = pinBin(bin_xi',stabXqs,stabXqs,'left');

    if sum(pSy.n >= cbin_min & ~isnan(pSy.qm)) >= nBin
        pSy.qm(pSy.n < cbin_min) = NaN;
        
        % stable_side_fit
        xref = linspace(0,1e-1);
        lfm  = fitlm(pSx.qm,pSy.qm,'linear','RobustOpts','on',...
                     'Intercept',itc_opt);        
        if Fitc
            incpt = lfm.Coefficients.Estimate(1);
            slp   = lfm.Coefficients.Estimate(2);
        else
            incpt = 0;
            slp   = lfm.Coefficients.Estimate(1);
        end
        [yref,yCI] = predict(lfm,xref');
        
        pfitStr = ['slope = ',num2str(round(slp,2),'% 4.2f'),...
                 ', intercept = ',num2str(round(incpt,4),'% 6.4f')];        
        fill(axm,[xref flip(xref)]',[yCI(:,1);flip(yCI(:,2))],...
                 rgb('coral'),'LineStyle','none','FaceAlpha',.1)
        pfit = plot(axm,xref,yref,'Color',rgb('coral'),'lineWidth',1.8);
        scatter(axm,pSx.qm,pSy.qm,200,'filled',...
                    'MarkerFaceColor',rgb('yellow'))
        errorbar(axm,pSx.qm,pSy.qm,pSy.qmL,pSy.qmU,'o','MarkerSize',13,...
                 'CapSize',10,'Color',rgb('ultramarine'),...
                 'linewidth',1.2,'linestyle','none')
    end
end

axm.YAxis.Exponent = -3;
axm.XAxis.Exponent = -3;
axm.Position = axm.Position + off_r*[1 1 -1 -1];
xlh = xlabel(axm,xltext,'fontsize',18);
ylh = ylabel(axm,yltext,'fontsize',18);
xlh.Units = 'Normalized';
ylh.Units = 'Normalized';
set(xlh,'Position',xlh.Position + [0 -off_r/2 0]);
set(ylh,'Position',ylh.Position + [-off_r/2 0 0]);

plot(axm,[-.1 .2],[-.1 .2],'--','color',.7*ones(1,3),'linewidth',2)
plot(axm,[-.1 .2],[0 0],':k'); plot(axm,[0 0],[-.1 .2],':k')
grid(axm,'on'); axis(axm,'square')
set(axm,'TickDir','out','Box','on')
xlim(axm,[-2e-3 2.2e-2]); ylim(axm,[-2e-3 2.2e-2])

legend(axm,pfit,pfitStr,'location','southeast','fontsize',16,...
       'AutoUpdate','off')

end

%% Subfunction for plotting horizontal pdf

function plot_xpdf(axPdf,pdfData,cStr)

histf = figure;
Shis  = histogram(pdfData,'Normalization','pdf','BinWidth',0.0005);
pdf   = Shis.Values;
bni   = Shis.BinEdges;
bn    = (bni(1:end-1)+bni(2:end))/2;

plot(axPdf,bn,pdf,'color',rgb(cStr),'linewidth',0.8)
hold(axPdf,'on')
fill(axPdf,[bn flip(bn)],[zeros(1,Shis.NumBins) flip(pdf)],...
     rgb(cStr),'LineStyle','none','FaceAlpha',.15)
set(axPdf,'XTick',[],'YTick',[],'Color','none','Box','off')

close(histf)
end

%% Subfunction for plotting vertical pdf

function plot_ypdf(axPdf,pdfData,cStr)

histf = figure;
Shis  = histogram(pdfData,'Normalization','pdf','BinWidth',0.0005);
pdf   = Shis.Values;
bni   = Shis.BinEdges;
bn    = (bni(1:end-1)+bni(2:end))/2;

plot(axPdf,pdf,bn,'color',rgb(cStr),'linewidth',0.8)
hold(axPdf,'on')
fill(axPdf,[zeros(1,Shis.NumBins) flip(pdf)],[bn flip(bn)],...
     rgb(cStr),'LineStyle','none','FaceAlpha',.15)
set(axPdf,'YTick',[],'Color','none','Box','off')

close(histf)
end

end