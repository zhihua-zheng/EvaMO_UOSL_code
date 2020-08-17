function [SKF,PROF,dprofT,dprofS,zSt,ziSt] = get_skf_prof(Lave,casename)

%% Parsing inputs

data_dir = './data/';
cLave = num2str(Lave);

switch casename
    
    case 'Papa'
        data_dir = [data_dir,'Papa/'];
        PROFname = fullfile(data_dir,['ocsp_prof_',cLave,'hrPMEL.mat']);
        FLUXname = fullfile(data_dir,['ocsp_flux_',cLave,'hrPMEL.mat']);
        WAVEname = fullfile(data_dir,['ocsp_wave_',cLave,'hr.mat']);
        lat      = 50.1;
        lon      = -144.9;
        ssp      = 1;
        
    case 'SPURSI'
        data_dir = [data_dir,'SPURSI/'];
        PROFname = fullfile(data_dir,['spursi_prof_',cLave,'hrBox.mat']);        
        FLUXname = fullfile(data_dir,['spursi_flux_',cLave,'hrUOP.mat']);
        WAVEname = fullfile(data_dir,['spursi_wave_',cLave,'hr.mat']);
        lat      = 24.5811;
        lon      = -38;
        ssp      = 0.865;
end

P = load(PROFname);
F = load(FLUXname);
W = load(WAVEname);

zSt     = W.zSt;
zSt_mid = (zSt(1:end-1) + zSt(2:end))/2;
ziSt    = union(zSt(1:end-1),zSt_mid);

dprofT = P.depth_t;
dprofS = P.depth_s;

%% Time parsing

datm_s = max([W.SD.datm(1)   F.SF.datm(1)   P.PROF.datm(1)  ]);
datm_e = min([W.SD.datm(end) F.SF.datm(end) P.PROF.datm(end)]);
idatm  = (datm_s:hours(Lave):datm_e)';
idahr  = datenum(idatm)*24; % in hours
ntm    = length(idatm);

InP = P.PROF.dahr >= idahr(1) & P.PROF.dahr <= idahr(ntm);
InM = F.SF.dahr   >= idahr(1) & F.SF.dahr   <= idahr(ntm);
InW = W.SD.dahr   >= idahr(1) & W.SD.dahr   <= idahr(ntm);

PROF = P.PROF(InP,:);
SF   = F.SF(InM,:);
SD   = W.SD(InW,:);

% find night time indices based on solar radiation, excluding transition
% periods (edges)

mNSW   = mode(SF.nsw); % night time nsw
allN   = find(SF.nsw < mNSW + .01);
d_allN = diff(allN); % index difference for night time
NRedge = [find(d_allN>1); length(allN)]; % right edge of night window
NLedge = [1;          find(d_allN>1)+1]; % left  edge of night window
Nedges = union(NRedge,NLedge); % all edges of night window
Ninter = setdiff(1:length(allN),Nedges); % all inter night

SKF.Idawn   = ismember(1:ntm,allN(NRedge))';
SKF.Idusk   = ismember(1:ntm,allN(NLedge))';
SKF.Inighte = ismember(1:ntm,allN(Nedges))';
SKF.Inighti = ismember(1:ntm,allN(Ninter))';
SKF.Idayi   = not(ismember(1:ntm,allN))';

% plot(SF.datm,SF.nsw);hold on;
% plot(SF.datm(SKF.Idawn),SF.nsw(SKF.Idawn),'.')
% plot(SF.datm(SKF.Inighte),SF.nsw(SKF.Inighte),'.');axis tight

% season separation
idvec       = datevec(idatm);
SKF.Imon    = idvec(:,2);
SKF.Iwinter = SKF.Imon == 12 | SKF.Imon == 1  | SKF.Imon == 2;
SKF.Ispring = SKF.Imon == 3  | SKF.Imon == 4  | SKF.Imon == 5;
SKF.Isummer = SKF.Imon == 6  | SKF.Imon == 7  | SKF.Imon == 8;
SKF.Iautumn = SKF.Imon == 9  | SKF.Imon == 10 | SKF.Imon == 11;

% season vector
SKF.Vseason              = zeros(ntm,1);
SKF.Vseason(SKF.Iautumn) = 1;
SKF.Vseason(SKF.Iwinter) = 2;
SKF.Vseason(SKF.Ispring) = 3;
SKF.Vseason(SKF.Isummer) = 4;

%% Constants

g     = 9.81;
rho0  = 1025;
cp    = 3985;
kappa = 0.4;

%% Seawater expansion/contraction coefficients

% ignore the depth-dependence of alpha and beta

% sea surface SA and CT before calibration
ssSA = gsw_SA_from_SP(SF.sss,ssp,lon,lat); % [g/kg]
ssCT = gsw_CT_from_t(ssSA,SF.sst,ssp);

[~,alpha,beta] = gsw_specvol_alpha_beta(ssSA,ssCT,ssp);

%% Kinematic fluxes due to turbulent processes

% net surface heat flux Qtur into the ocean
Qtur = SF.hlb + SF.hsb + SF.nlw; % [W/m^2]

% bulk surface heat flux
SKF.Qo = SF.hlb + SF.hsb + SF.nlw + SF.nsw; % [W/m^2]

% surface freshwater flux into the ocean
Ftur = (SF.rain - SF.evap)/1000/3600; % [m/s]

% surface kinematic fluxes
SKF.w_theta_0 = -Qtur/cp/rho0; % <w't'>_0 [C*m/s]
SKF.w_s_0     =  Ftur.*ssSA;   % <w's'>_0 [(g/kg)*(m/s)]
SKF.wb0T      =  g*(alpha .* SKF.w_theta_0);
SKF.wb0S      = -g*(beta  .* SKF.w_s_0);
SKF.w_b_0     =  SKF.wb0T + SKF.wb0S; % <w'b'>_0 [m^2/s^3]

SKF.Ustar = sqrt(SF.tau/rho0); % [m/s]

% bulk surface buoyancy flux, positive for buoycancy gain
SKF.Bo = -g*(alpha .* (-SKF.Qo/cp/rho0) - beta .* SKF.w_s_0); % [m^2/s^3]

SKF.u10   = SF.u10;
SKF.nsw   = SF.nsw;
SKF.alpha = alpha;
SKF.beta  = beta;
SKF.rain  = SF.rain;
SKF.evap  = SF.evap;

%% Stokes drift velocity

SKF.lstd  = SKF.u10.^2/g; % e-folding scale for waves
SKF.SD    = SD;
SKF.USt_0 = sqrt(SD.uSt(:,end).^2 + SD.vSt(:,end).^2);

%% Stokes drift shear projected to wind direction

Ustar_kz  = SKF.Ustar / kappa ./ abs(ziSt)';
windUnitx = SF.tau_x ./ SF.tau;
windUnity = SF.tau_y ./ SF.tau;

% downwind component & crosswind component of Stokes drift shear
SKF.dUStDw_dz = windUnitx .* SD.duSt_dz + windUnity .* SD.dvSt_dz;
SKF.dUStCw_dz = windUnitx .* SD.dvSt_dz - windUnity .* SD.duSt_dz;

SKF.UStDw = windUnitx .* SD.uSt + windUnity .* SD.vSt;
SKF.UStCw = windUnitx .* SD.vSt - windUnity .* SD.uSt;

SKF.eta_x = -SKF.dUStDw_dz ./ Ustar_kz;
SKF.eta_y = -SKF.dUStCw_dz ./ Ustar_kz;

% direction of wind and Stokes dirft, counter-clockwise from East
SKF.wind_deg = atan2d(windUnity,windUnitx);
SKF.wave_deg = atan2d(SD.vSt,SD.uSt);

%% roughness length

zb = 0.6; % z0/Hs
SKF.z0 = zb*SD.Hs;

end