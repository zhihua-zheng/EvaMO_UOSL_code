function [gradT,PTfit,itype] = get_fit_gdT(depth,iD,PTprof,BLD,nFit,ishow)

%% Pre-setting

lnz = log(depth);
[~,ntm] = size(PTprof);

if iD < 3
    rfit = [0.1 0.3];
else
    rfit = [0.2 0.3];
end

Dfit = BLD*rfit;
SLD  = BLD*0.2;
nL   = find(depth < max(SLD),1,'last');

gradT = nan(nL,ntm);
PTfit = nan(nL,ntm);
itype = zeros(1,ntm); % id for fitting type
% itype = 0, insufficient data points
% itype = 1, fit is not robust
% itype = 2, successful fit

%% Second-order polynomial fit of profile in depth range [0, 0.3H]

for j = 1:ntm
    
    % number of levels available for least-square fitting
    [~,mlim] = min(abs(depth - Dfit(j,:))); % 
    m  = (mlim(1):mlim(2));
    nr = length(m);
    
    % number of measurements within surface layer
    iSLDj   = sum(depth <= SLD(j));
    nPTinSL = sum(~isnan(PTprof(1:mlim(2),j)));
    
    %% loop for regressions with varying points

    if nPTinSL >= nFit
        
        p   = cell(1,nr);
        gof = cell(1,nr);
        pTex     = nan(3,nr);
        rms_pTex = nan(1,nr);
        Imono    = zeros(1,nr);
        I2mono   = zeros(1,nr);
        Isense   = zeros(1,nr);
        Ismall   = zeros(1,nr);
        
        for ir = 1:nr
            Ifit = ~isnan(PTprof(1:m(ir),j));

            if sum(Ifit) >= 3
                [p{ir},gof{ir}] = fit(lnz(Ifit),PTprof(Ifit,j),'poly2');

                % examine extrapolation accuracy
                pTex(:,ir) = p{ir}(lnz(m(ir):m(ir)+2));
                rms_pTex(ir) = rms(pTex(:,ir) - PTprof(m(ir):m(ir)+2,j));

                % examine the gradient change of the fitting function
                dz   = 0.1;
                zref = log(depth(1):dz:depth(m(ir)));
                dpT  = diff(p{ir}(zref),1);
                d2pT = diff(p{ir}(zref),2);
                ieG  = abs(dpT)  >= 1e-7; % ignore the sign of small dT/dz
                i2eG = abs(d2pT) >= 1e-8; % ignore the sign of small d2T/dz2
                sign_dpT  = sign(dpT(ieG));
                sign_d2pT = sign(d2pT(i2eG));
                Imono(ir)  = length(unique(sign_dpT))  == 1;
                I2mono(ir) = length(unique(sign_d2pT)) == 1;
                Isense(ir) = gof{ir}.adjrsquare > 0.5;
                Ismall(ir) = gof{ir}.rmse < 2e-3;
            end
        end
        
        % good fit: monotonic in second-order, large r2, small rmse
        igood = find(Imono & I2mono & Isense & Ismall);
        if isempty(igood)
            itype(j) = 1;
            continue
        end
        
        % best fit: miminum rmse
        [~,ibest] = min(cellfun(@(A) A.rmse,gof(igood)));
        bestR = igood(ibest);
        M     = m(bestR);
        P     = p{bestR};
        Gof   = gof{bestR};
        
        % polynomial coefficients
        b = [p{bestR}.p1 p{bestR}.p2 p{bestR}.p3];
        
        %% compute vertical gradient & fitted temperature
                    
        % close to SL base, gradient is small, susceptible to error
        nGrad   = min([M-1 iSLDj]);
        Ieva    = depth(1:nGrad) > depth(1);
        d_eva   = depth(Ieva);
        lnz_eva = log(d_eva);
        gradT(Ieva,j) = -(2*b(1)*lnz_eva + b(2)) ./ d_eva;

        % fitted profile in surface layer
        PTfit(1:iSLDj,j) = P(lnz(1:iSLDj));
        itype(j) = 2;

        % visualize polynomial regression
        if ishow
            show_polyfit_lnz;

        % Figure 3 is from the SPURS-I dataset with j = 3427
        elseif iD == 2 && j == 3427
            show_polyfit_lnz;
        end                
    end
end

%% Subfunction for visualization of regression

function show_polyfit_lnz()

    zref = linspace(lnz(1),lnz(M));
    pT = P(zref);
    off_r = 0.02;

    figure('position',[10 30 695 620])
    [ax,~] = tight_subplot(1,2,[.02 .05],[.1 .05],[.1 .1]);

    plot(ax(1),PTprof(:,j),-lnz,'o--','markersize',6); hold(ax(1),'on');
    plot(ax(1),pT,-zref,'linewidth',2)
    plot(ax(1),pTex(:,bestR),-lnz(M:M+2),'-.r','linewidth',2)
    plot(ax(1),[PTprof(1,j) max(PTprof(1:M+2,j))],-log(SLD(j))*[1 1],'--k')
    scatter(ax(1),PTprof(Ieva,j),-lnz_eva,35,rgb('goldenrod'),'filled')
    grid(ax(1),'on')
    set(ax(1),'TickDir','out','ylim',[-lnz(M+2) -lnz(1)+0.1])
    text(ax(1),0.98,0.98,'(a)','Units','Normalized','FontSize',22,...
        'HorizontalAlignment','right','VerticalAlignment','top')
    text(ax(1),0.35,0.8,['RMSE: ',num2str(round(Gof.rmse,4)),' ',...
        char(0176),'C'],'fontsize',16,'Units','Normalized')
    ax(1).Position = ax(1).Position + off_r*[1 1 -1 -1];
    xlh = xlabel(ax(1),['Temperature [',char(0176),'C]'],'fontsize',18);
    ylh = ylabel(ax(1),'-ln(|z|) [m]','fontsize',18);
    xlh.Units = 'Normalized';
    ylh.Units = 'Normalized';
    set(xlh,'Position',xlh.Position + [0 -off_r 0]);
    set(ylh,'Position',ylh.Position + [-off_r 0 0]);

    plot(ax(2),PTprof(:,j),-depth,'o--','markersize',6); hold(ax(2),'on');
    plot(ax(2),pT,-exp(zref),'linewidth',2)
    plot(ax(2),pTex(:,bestR),-depth(M:M+2),'-.r','linewidth',2)
    plot(ax(2),[PTprof(1,j) max(PTprof(1:M+2,j))],-SLD(j)*[1 1],'--k')
    scatter(ax(2),PTprof(Ieva,j),-depth(Ieva),35,rgb('goldenrod'),'filled')
    grid(ax(2),'on')
    set(ax(2),'TickDir','out','YAxisLocation','right','ylim',[-depth(M+2) 0])
    text(ax(2),0.98,0.98,'(b)','Units','Normalized','FontSize',22,...
        'HorizontalAlignment','right','VerticalAlignment','top')
    text(ax(2),0.4,0.7,['r_*^2: ',num2str(round(Gof.adjrsquare,2))],...
        'fontsize',16,'Units','Normalized')
    ax(2).Position = ax(2).Position + off_r*[0 1 -1 -1];
    xlh = xlabel(ax(2),['Temperature [',char(0176),'C]'],'fontsize',18);
    ylh = ylabel(ax(2),'z [m]','fontsize',18);
    xlh.Units = 'Normalized';
    ylh.Units = 'Normalized';
    set(xlh,'Position',xlh.Position + [0 -off_r 0]);
    set(ylh,'Position',ylh.Position + [off_r 0 0]);
end

end