
%% Load processed datasets

root_dir   = '~/GDrive/UW/Research/Data/';
ocsp_dir   = [root_dir,'Papa/'];
spursi_dir = [root_dir,'SPURSI/'];

pzPA = load([ocsp_dir,  'ocsp_pzData.mat']);
pzSP = load([spursi_dir,'spursi_pzData.mat']);

skfPA = load([ocsp_dir,  'ocsp_skfData.mat']);
skfSP = load([spursi_dir,'spursi_skfData.mat']);
PAqs  = skfPA.Islc;
SPqs  = skfSP.Islc;

dTPA = load([ocsp_dir,  'ocsp_dTData.mat']);
dTSP = load([spursi_dir,'spursi_dTData.mat']);

rpPA = load([ocsp_dir,  'ocsp_rpData.mat']);
rpSP = load([spursi_dir,'spursi_rpData.mat']);

%% Figure 4: temperature comparison between observations and MOST

Poff = 0.055;
dTf = figure('position',[0 0 800 800]);
[dTax,dTpos] = tight_subplot(2,2,[.01 .02],[.06 .05],[.06 .08]);
dTax(1).Position = [dTpos{1}(1:2) dTpos{1}(3)+Poff dTpos{1}(4)];
dTax(2).Position = [dTpos{2}(1)+Poff dTpos{2}(2) dTpos{2}(3)-Poff dTpos{2}(4)];
dTax(3).Position = [dTpos{3}(1:2) dTpos{3}(3)+Poff dTpos{3}(4)];
dTax(4).Position = [dTpos{4}(1)+Poff dTpos{4}(2) dTpos{4}(3)-Poff dTpos{4}(4)];

dT_scat(dTPA,PAqs,1,dTax(1),1,1);
dT_scat(dTSP,SPqs,2,dTax(3),1,1);
plot_profcom(rpPA,1,dTax(2));
plot_profcom(rpSP,2,dTax(4));

%% Figure 5: boxplot for observed phi_h in zeta space

phi_obs = {pzPA.phi_qs(:)  pzSP.phi_qs(:)};
zet_obs = {pzPA.zeta_qs(:) pzSP.zeta_qs(:)};
PZ_scat_obs(phi_obs,zet_obs);

%% Figure 6: predictions from the 'SE' surface wave breaking model

Poff = 0.08;
WBf = figure('position',[0 0 935 530]);
[WBax,WBpos] = tight_subplot(1,2,[.02 .05],[.09 .1],[.05 .09]);
WBax(1).Position = [WBpos{1}(1:2) WBpos{1}(3)+Poff WBpos{1}(4)];
WBax(2).Position = [WBpos{2}(1)+Poff WBpos{2}(2) WBpos{2}(3)-Poff WBpos{2}(4)];

r2o = -0.2; % z_0/L
plot_SFWB_phi2d(WBax(1),r2o,pzPA,pzSP)
q_in_SFWB(r2o,[],WBf,WBax(2))

%% Figure 7: predictions from the SE Langmuir turbulence model

figure('position',[0 0 1015 500]);
[LTax,~] = tight_subplot(1,2,[.03 .03],[.1 .1],[.06 .08]);
plot_H15se_phi2d(LTax,'H2015',pzPA,pzSP);

%% Figure 8: comparison of phi_h from Kansas, observations and models

PZf = figure('position',[0 0 1050 510]);
[PZax,~] = tight_subplot(1,2,[.01 .01],[.12 .03],[.08 .03]);

PZ_scat(pzPA,1,PZf,PZax(1),1);
PZ_scat(pzSP,2,PZf,PZax(2),1);

%% Figure 9: phi_h & l^* from the SE Langmuir turbulence model including length scale Eq.

% Find optimal E6 value for SMCLT with varying l^*

% dT_obs_SMCLT;
% figure;
% plot(E6_list,slp,'-o'); grid on
% optimal E6 is 2.4, with slp = 1
save([spursi_dir,'spursi_slpE6.mat'],'E6_list','slp')

Lstarf = figure('position',[0 0 820 800]);
[lsrAx,~] = tight_subplot(2,2,[.07 .02],[.07 .03],[.08 .03]);

phi_lstar_E6(pzPA,1,oE6,lsrAx,1);
phi_lstar_E6(pzSP,2,oE6,lsrAx,1);

%% Forcing distribution at two sites

xString = 'H_s [m]';
plot_pdf(skfPA.SKF.SD.Hs(PAqs),skfSP.SKF.SD.Hs(SPqs),22,xString)
xlim([0 10])

xString = 'u_* [m/s]';
plot_pdf(skfPA.SKF.Ustar(PAqs),skfSP.SKF.Ustar(SPqs),22,xString)
xlim([0 0.04])

xString = 'B_0 [m^2/s^3]';
plot_pdf(skfPA.SKF.Bo(PAqs),skfSP.SKF.Bo(SPqs),42,xString)
xlim([-6e-7 4e-7])

xString = '\eta^x';
plot_pdf(pzPA.etaX_qs(:),pzSP.etaX_qs(:),100,xString)
xlim([-1.3 .3])

xString = '\eta^y';
plot_pdf(pzPA.etaY_qs(:),pzSP.etaY_qs(:),120,xString)
xlim([-1 1])
