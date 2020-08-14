function [rPTprof,rPTmost,zml,Iprof] = spursi_rp(PTprof,depth_t,...
                                       BLD,dataName,SKF,Islc)

SLD = BLD/5;

% find the mean temperature in surface layer
mPTprof = get_mT_inLayer(PTprof,depth_t,SLD,@nanmean);
rPTprof = PTprof - mPTprof; % relative to mPTprof
zml     = -depth_t ./ BLD';

% temperature profile from MOST integration
zdum = -[.75 (1:.1:20) (21:101) (105:5:200) 200.5 (250:50:400)]' - 0.115;
[MOL,zeta,FS] = MOSTpar_from_flux(zdum,dataName,SKF);
PTmost  = get_MOST_prof(-depth_t,zdum,PTprof(1,:),FS.Tstar,zeta,'Kansas');
mPTmost = get_mT_inLayer(PTmost,depth_t,SLD,@nanmean);
rPTmost = PTmost - mPTmost; % relative to mPTmost

%% Indices for unstable side

IspanSL = max(zml) > -0.02 & min(zml) < -0.2;
Ineg    = sum(MOL>0) == 0;
Iprof   = Islc & Ineg & IspanSL;

end