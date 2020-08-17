%% prep_spursi
%
% Process data from the SPURS-I site
%
% Zhihua Zheng, University of Washington, July 18 2019
% =========================================================================

%% Constants

clear
Lave = 1; % the length of averaging interval [hr]
kappa = 0.4;

%% Loading

iD = 2;
scat_mode = 1;
dataName = 'SPURSI';
data_dir = './data/';
spursi_dir = [data_dir,'SPURSI/'];

[SKF,PROF,depth_t,depth_s,zSt,ziSt] = get_skf_prof(Lave,dataName);

%% observed profiles

idatm = PROF.datm;
idahr = PROF.dahr;
itime = datenum(idatm);
iday  = itime - itime(1);
idatm_loc = idatm - hours(timezone(-38));

% calibrate salinity
cPROF = calib_spursi_S(PROF,depth_t,depth_s,SKF);

PTprof  = cPROF.PTprof';
SAprof  = cPROF.SAprof';
SAproff = cPROF.SAproff';
PDprof  = cPROF.PDprof';
Bprof   = cPROF.Bprof';
NSQprof = cPROF.NSQprof';

%% Richadrson number & BLD

nziSt    = length(ziSt);
ziSt_up  = flip(ziSt);
[nz,ntm] = size(PTprof);
zmid     = -(depth_t(2:end) + depth_t(1:end-1))/2;

% gradient & bulk Richardson number, using idealized shear in log layer
SSQprof = (repmat(SKF.Ustar',nz-1,1) ./ repmat(zmid,1,ntm) / kappa).^2;
Rigprof = NSQprof ./ SSQprof;
Ribprof = get_Rib(Bprof,SKF.Ustar',-depth_t);

% MLD = get_mld(flipud(PDprof),-flipud(depth_t),1,-1);
% Rib_c_list = (0.25:0.05:8);
% BLD_rmse = nan(length(Rib_c_list),1);
% for i = 1:length(Rib_c_list)
%     BLD = get_mld(flipud(Ribprof),-flipud(depth_t(2:end)),3,Rib_c_list(i));
%     BLD_rmse(i) = nanrms(BLD-1.2*MLD);
% end
% figure('position',[0 0 600 300])
% plot(Rib_c_list,BLD_rmse,'-o'); grid on

% 2.30 for 1.0MLD
BLD = get_mld(flipud(Ribprof),-flipud(depth_t(2:end)),3,2.30);
SLD = BLD*0.2;
nL  = find(depth_t < max(SLD),1,'last');
ziSt_SL = ziSt_up(1:find(ziSt_up == -depth_t(nL)));

% position of SL sensor depths in ziSt_up
[~,Pct] = ismember(-depth_t(1:nL),ziSt_up);

% indices of SL sensor depths in ziSt_up
Ict = ismember(ziSt_up,-depth_t(1:nL));

% indices of SL sensor depths in ziSt_SL
ISL_ct = ismember(ziSt_SL,-depth_t(1:nL));

%% Polynomial fit of profiles in ln(|z|)

[gdT_fit,PTfit] = get_fit_gdT(depth_t,iD,PTprof,BLD,3,0);
save([spursi_dir,'spursi_fitData.mat'],'gdT_fit','PTfit')

%% Surface proximity function

Ls  = get_Ls(ziSt,BLD,SKF.dUStDw_dz,kappa);
fzS = 1 + tanh(ziSt/4./Ls);
SKF.Ls = Ls';

%% Normalized Stokes drift shear and fzs in Harcourt (2015)

% transposition
etaX = flipud(SKF.eta_x');
etaY = flipud(SKF.eta_y');
fzS  = flipud(fzS);

%% Partial convection and quasi-steady state

% heat flux can change from negative at surface to positive at depth,
% due to the penetrative shortwave radiation

[MOL_SL,zeta_SL,FS_SL] = MOSTpar_from_flux(ziSt_SL,dataName,SKF);
IpartConvec = FS_SL.Tstar(1,:) < 0 & FS_SL.Tstar(end,:) > 0;

% buoyancy flux at boundary layer depth
[~,~,FS_BLD] = MOSTpar_from_flux(-BLD,dataName,SKF);
BfH = FS_BLD.Bf;
get_qs_time;

Iqs  = (stage == 1)'; % quasi-steady state
Islc = and(Iqs,~IpartConvec);
nSlc = sum(Islc);
save([spursi_dir,'spursi_skfData.mat'],'SKF','Islc')

%% Vertical temperature differences from observations & MOST

% vertical temperature differences
dT_obs = nan((nL*(nL-1))/2,ntm);
dT_MOi = nan(size(dT_obs));
dT_fit = nan(size(dT_obs));

ic  = 0;
fop = 'Kansas';

for iu = 1:nL-1      % index for upper level
    for jl = iu+1:nL % index for lower level

ic = ic + 1;
zz = -depth_t([iu,jl]);

dT_obs(ic,:) = PTprof(iu,:) - PTprof(jl,:);
dT_fit(ic,:) = PTfit(iu,:)  - PTfit(jl,:);

[MOL,~,FS] = MOSTpar_from_flux(zz,dataName,SKF);

% integration within two levels to estimate temperature difference
dT_MOi(ic,:) = get_MOST_Delta(FS.Tstar,MOL,zz,fop);

% only use data in the surface layer
IinSL = -SLD' < zz(2);
dT_obs(ic,~IinSL) = NaN;
dT_MOi(ic,~IinSL) = NaN;
    end
end

save([spursi_dir,'spursi_dTData.mat'],'dT_MOi','dT_fit','dT_obs')

%% Temperature profiles for comparsion

[rObs,rMOi,zbl,Iprof] = spursi_rp(PTprof,depth_t,BLD,dataName,SKF,Islc);
save([spursi_dir,'spursi_rpData.mat'],'rObs','rMOi','zbl','Iprof')

%% Compute dimensionless temperature gradients

gdT_log = FS_SL.Tstar(ISL_ct,:) ./ (-ziSt_SL(ISL_ct));
phi_fit = gdT_fit ./ gdT_log;
zet_ct  = depth_t(1:nL) ./ MOL_SL(ISL_ct,:);
xi_ct   = depth_t(1:nL) ./ SKF.z0';

etaXct = etaX(Ict,:);
etaYct = etaY(Ict,:);
fzSct  =  fzS(Ict,:);

% Save PZ data
phi_qs  = phi_fit(:,Islc);
zeta_qs = zet_ct(:,Islc);
etaX_qs = etaXct(:,Islc);
etaY_qs = etaYct(:,Islc);
fzS_qs  = fzSct(:,Islc);
xi_qs   = xi_ct(:,Islc);

save([spursi_dir,'spursi_pzData.mat'],...
     'phi_qs','zeta_qs','etaX_qs','etaY_qs','fzS_qs','xi_qs')
clear