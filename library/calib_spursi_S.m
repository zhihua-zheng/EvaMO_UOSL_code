function cPROF = calib_spursi_S(iPROF,depth_t,depth_s,SKF)

%% Loading

idatm   = iPROF.datm;
itime   = datenum(idatm);
iday    = itime - itime(1);
iPTprof = iPROF.PTprof';
iSAprof = iPROF.SAprof';
iPDprof = iPROF.PDprof';
tprof   = iPROF.tprof';
w_s_0   = SKF.w_s_0;
w_t_0   = SKF.w_theta_0; % night time, no solar radiation

Idawn   = find(SKF.Idawn);
Idusk   = find(SKF.Idusk);
Inighti = find(SKF.Inighti);

%% Exclude the hours affected by rain events for each night

[nz,ntm] = size(iPTprof);
nDawn    = length(Idawn);

nNi = sum(SKF.Inighti);
for j = 1:nNi   
    Ind = Inighti(j);
    if sum( SKF.rain(Ind-3:Ind)>0 ) > 0 % rained in the past 3 hours
        Inighti(j) = NaN;
    end
end

Inighti(isnan(Inighti)) = [];
nNighti = length(Inighti);

%% Depth of the convective layer from temperature profile

PTnighti = iPTprof(:,Inighti);

% 0.006 C criteria
dPT    = 0.006;
ref_d  = 3;
ref_PT = nanmean(PTnighti(depth_t<ref_d,:));

% convective layer depth, 0.006C warmer/colder than surface
where_utd = PTnighti - ref_PT;
Jutd = nan(nNighti,1);

for n = 1:nNighti   
    col = abs(where_utd(:,n));
    if col(1) <= dPT
        Jutd(n) = find(col > dPT,1,'first') - 1;
    end
end

%% Find the deepest homogeneous temperature profile for each night

% indices for selected convective temperature profile
Ict = nan(size(Idawn));
Jct = nan(size(Idawn));

for n = 1:nDawn
    
    % indices for each night
    IIpool = Inighti > Idusk(n) & Inighti < Idawn(n);
    
    if sum(IIpool) > 0        
        PTpool = PTnighti(:,IIpool); % temp. profile in that night
        Ione_night = Inighti(IIpool); % indices for that night in ntm
        Jutd_pool  = Jutd(IIpool); % utd indices for that night
        [mJutd,IIIdeep] = max(Jutd_pool); % deepest profile
        Jct(n) = mJutd; % utd indices for choosen profile

        if length(IIIdeep) > 1 % minimum variance
            [~,wellmix] = min(var(PTpool(1:mJutd,IIIdeep)));
            Ict(n) = Ione_night(IIIdeep(wellmix));
        else
            Ict(n) = Ione_night(IIIdeep);
        end
    end
end

ct_datm = NaT(size(Ict));
ind = ~isnan(Ict);
ct_datm(ind) = idatm(Ict(ind));
ct_time = datenum(ct_datm);
ct_day  = ct_time - itime(1);

% plot(ct_datm(ind),SKF.rain(Ict(ind)),'.-')
% ylabel('rain [mm/hr]'); axis tight

%% Shrink the salinity variation in convective layer

SA_ict = nan(nz,nDawn);
PT_ict = nan(nz,nDawn);
SA_off_ict = nan(nz,nDawn);

SA_ict(:,ind) = iSAprof(:,Ict(ind)); % profile at Ict to be corrected
PT_ict(:,ind) = iPTprof(:,Ict(ind));
cSA_ict = SA_ict; % profile after correction

% The maximum of salinity variation (relative to the mean) within the
% convective layer is set to scale with that of temperature.
% The ratio is determined by the surface kinematic fluxes [see Eq. (B1)]

dPT_max = nan(nDawn,1);
dSA_max = nan(nDawn,1);
scale_fac = nan(nDawn,1);

for n = 1:nDawn
    
    jutd = Jct(n);
    
    if ~isnan(jutd)
        dPT_convec = PT_ict(1:jutd,n) - nanmean(PT_ict(1:jutd,n));
        dSA_convec = SA_ict(1:jutd,n) - nanmean(SA_ict(1:jutd,n));

        dSA_max(n) = max(abs(dSA_convec));
        dPT_max(n) = max(abs(dPT_convec));

        dSA_scale = dPT_max(n) * abs(w_s_0(Ict(n)) / w_t_0(Ict(n)));
        scale_fac(n) = dSA_scale/dSA_max(n);

        % only adjust if dSA_max is larger than dSA_scale
        if scale_fac(n) < 1
            dSA_shrink = scale_fac(n)*dSA_convec;
            cSA_ict(1:jutd,n) = nanmean(SA_ict(1:jutd,n)) + dSA_shrink;
            % value of adjustment
            SA_off_ict(1:jutd,n) = cSA_ict(1:jutd,n) - SA_ict(1:jutd,n);
        end
    end
end

%% Extract the slow-varying part of the offset time series

% 29 points median filter
SA_off_fil = medfilt1(SA_off_ict,29,[],2,'omitnan','truncate');
SA_off = nan(nz,ntm);

for n = [1 3 6 8:14]
    Icor = and(~isnan(SA_off_fil(n,:)),~isnat(ct_datm)');
    SA_off(n,:) = interp1(ct_datm(Icor),SA_off_fil(n,Icor),...
                          idatm,'linear','extrap');
end

for n = [15:16 18]
    SA_off_fil(n,:) = medfilt1(SA_off_ict(n,:),61,[],2,'omitnan','truncate');
    Icor = and(~isnan(SA_off_fil(n,:)),~isnat(ct_datm)');
    SA_off(n,:) = interp1(ct_datm(Icor),SA_off_fil(n,Icor),...
                          idatm,'linear','extrap');
end

Joff = find(sum(~isnan(SA_off_ict),2)>0,1,'last');
for n = [17 19:Joff]
    Icor = ~isnan(SA_off_ict(n,:));
    lm_SA_off = fitlm(ct_day(Icor),SA_off_ict(n,Icor),'Intercept',false);
    SA_off(n,:) = predict(lm_SA_off,iday);
end

SA_offu = SA_off + .03;
SA_offl = SA_off - .03;

%% Drift curves part I

dateStr  = {'01/11/2012','01/02/2013','01/05/2013','01/08/2013'};
ticklocs = datetime(dateStr,'InputFormat','dd/MM/yyyy');
figure('position',[0 0 820 700])
ic = 0;
[hax1,~] = tight_subplot(5,2,[.03 .02],[.1 .05],[.1 .05]);

for n = [1 3 6 8:14]
    
    ic = ic + 1;
    axes(hax1(ic))
    xref = [idatm' flip(idatm)'];
    yref = [SA_offu(n,:) flip(SA_offl(n,:))];
    fill(xref,yref,rgb('orange'),'LineStyle','none','FaceAlpha',.2); hold on
    h1 = plot(ct_datm,SA_off_ict(n,:),'o','MarkerSize',3,'color',rgb('water blue'));
    setMarkerColor(h1,'water blue',.2);
    plot(idatm,SA_off(n,:),'linewidth',2,'color',rgb('orange red'))
    plot([idatm(1)-days(10) idatm(end)+days(10)],[0 0],'--k','linewidth',.5)
    grid on; ylim([-.12 .12]); xtickformat('MMM-yyyy')
    set(gca,'fontsize',14,'XMinorTick','on','YMinorTick','on',...
            'TickDir','out','XLimSpec','Tight','XTick',ticklocs)
    if ic < 9
       set(gca,'XtickLabel',[])
    end
    if mod(ic,2) == 0
       set(gca,'YtickLabel',[])
    end
    text(ct_datm(7),.07,['sensor at ',num2str(round(depth_t(n),2)),' m'],...
         'fontsize',16)
end

[~,hsupY1] = suplabel('\Delta [g kg^{-1}]','y');
set(hsupY1,'fontsize',18,'Interpreter','tex','Units','Normalized')
set(hsupY1,'Position',hsupY1.Position + [0.02 0 0]);

%% Drift curves part II

% figure('position',[0 0 820 700])
% ic = 0;
% [hax2,~] = tight_subplot(5,2,[.03 .02],[.1 .05],[.1 .05]);
% 
% for n = (15:24)
%     ic = ic + 1;
%     axes(hax2(ic))
%     h1 = plot(ct_datm,SA_off_ict(n,:),'o','MarkerSize',3);hold on
%     setMarkerColor(h1,'azure',.15);
%     plot(idatm,SA_off(n,:),'linewidth',1.8,'color',rgb('orange red'))
%     plot([idatm(1)-days(10) idatm(end)+days(10)],[0 0],'--k','linewidth',.5)
%     grid on; ylim([-.1 .1]); xtickformat('MMM-yyyy')
%     set(gca,'fontsize',14,'XMinorTick','on','YMinorTick','on',...
%             'TickDir','out','XLimSpec','Tight','XTick',ticklocs)
%     if ic < 9
%        set(gca,'XtickLabel',[])
%     end
%     if mod(ic,2) == 0
%        set(gca,'YtickLabel',[])
%     end
%     text(ct_datm(7),.07,['sensor at ',num2str(round(depth_t(n),2)),' m'],...
%          'fontsize',16)
% end
% [~,hsupY2] = suplabel('\Delta [g/kg]','y');
% set(hsupY2,'fontsize',18,'Interpreter','tex','Units','Normalized')
% set(hsupY2,'Position',hsupY2.Position + [0.02 0 0]);

%% Adjust the salinity profiles with drift curves

SA_off(isnan(SA_off)) = 0;
Icfail = ~isnan(SA_off);
SAprof = iSAprof + SA_off;

% cmap_div = flipud( cbrewer('div','RdBu',21) );
% figure('position',[0 0 620 350])
% ax_pc = pcolor(idatm,-depth_t,SA_off);
% colormap(cmap_div); cbar = colorbar('eastoutside'); caxis([-.2 .2])
% set(ax_pc,'EdgeColor','none')
% grid on; box on; xtickformat('MMM-yyyy'); ylim([-100 0])
% ylabel('z [m]','fontsize',16)
% set(gca,'fontsize',11,'XMinorTick','on','YMinorTick','on','TickDir','out')
% ylabel(cbar,'$\overline{S}_{cor}$ - $\overline{S}_{ori}$ [$g/kg$]',...
%        'fontsize',16,'Interpreter','latex');

%% Figure: resccale factor

% figure('position',[0 0 620 350])
% plot(ct_datm,scale_fac,'linewidth',1.2)
% grid on; box on; xtickformat('MMM-yyyy'); ylim([0 0.4])
% set(gca,'fontsize',11,'XMinorTick','on','YMinorTick','on','TickDir','out')
% ylh = ylabel(['$\big| \Delta \overline{S} \big|_{m} \Big/ ',...
%        '\big| \Delta \overline{S} \big|_{m, \: obs}$'],...
%        'fontsize',16,'Interpreter','latex');
% rsf_pos = get(gca,'Position');
% set(gca,'Position',rsf_pos + 0.03*[1 0 -1 0]);
% ylh.Units = 'Normalized';
% set(ylh,'Position',ylh.Position + [-0.03 0 0]); 

%% Recompute potential temperature with corrected salinity

SAproff = vert_fill(SAprof,depth_s,idatm);
PTprof  = gsw_pt0_from_t(SAproff,tprof,depth_t);

%% Recompute potential density, buoyancy

g = 9.81;
rho0 = 1025;

CTprof = gsw_CT_from_t(SAproff,tprof,depth_t);
PDprof = gsw_sigma0(SAproff,CTprof);
Bprof  = -g*(PDprof - rho0)/rho0;

% buoyancy frequency squared [1/s^2]
NSQprof = center_diff(Bprof,-depth_t,1,'mid');

%% Figure: example demonstrating the changes after correction

demo_spursi_cprof;

%% Timetable

PTprof  = PTprof';
SAprof  = SAprof';
SAproff = SAproff';
PDprof  = PDprof';
Bprof   = Bprof';
NSQprof = NSQprof';
Icfail  = Icfail';

iPROF.PTprof  = PTprof;
iPROF.SAprof  = SAprof;
iPROF.SAproff = SAproff;
iPROF.PDprof  = PDprof;
iPROF.Bprof   = Bprof;
iPROF.NSQprof = NSQprof;

cPROF = addvars(iPROF,Icfail);

end