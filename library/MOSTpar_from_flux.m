function [MOL,zeta,FS] = MOSTpar_from_flux(zz,casename,SKF)
%
% MOSTpar_from_flux
%==========================================================================
%
% USAGE:
%  [MOL,zeta,FS] = MOSTpar_from_flux(zz,casename,SKF)
%
% DESCRIPTION:
%  Compute parameters associated with Monin-Obukhov similarity theory from 
%  hourly time series of surface kinematic fluxes.
%
% INPUT:
%
%  idatm - 1-D column vector of MATLAB datetime for the inquired period
%  zz - vertical coordinates for the choosen 2 levels [-, m]
%  casename - string to indicate the parameters for which cases is inquired
%         
% OUTPUT:
%
%  MOL - Monin-Obukhov length for evaluated depth [m]
%  zeta - Monin-Obukhov stability parameter
%  FS - struct contains turbulence scale for velocity (Ustar), temperature 
%       (Tstar), salinity (Sstar), buoyancy (Bstar), buoyancy forcing Bf
%
% AUTHOR:
%  June 29 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Parsing inputs

switch casename
    
    case 'Papa'
        waterType = 4;
        
    case 'SPURSI'
        waterType = 3;
end

w_b_0     = SKF.w_b_0';
w_s_0     = SKF.w_s_0';
w_theta_0 = SKF.w_theta_0';
alpha     = SKF.alpha';
NSW       = SKF.nsw;

%% Constants

ntm = length(w_b_0);

MOconsts_name = '~/GDrive/UW/Research/Data/Misce/MO_consts.mat';
load(MOconsts_name,'g','rho0','cp','kappa');

%% T-S dependent coefficients

% ignore the depth-dependence as T-S variation with depth is small

% use sea surface SA and CT before calibration
% ssSA = gsw_SA_from_SP(SKF.sss,ssp,lon,lat); % [g/kg]
% ssCT = gsw_CT_from_t(ssSA,SKF.sst,ssp);
% 
% % isobaric heat capacity of seawater [J/kg/C]
% % cp    = gsw_cp_t_exact(ssSA,SKF.sst,ssp);
% 
% alpha = gsw_alpha(ssSA,ssCT,ssp);
% % beta  = gsw_beta(ssSA,ssCT,ssp);

%% Fluxes due to shortwave radiation

band_SR = 2; % 2 band model
ndum    = 101;

if length(zz) == 2

    z1   = zz(1); % upper level
    z2   = zz(2); % lower level
    zdum = linspace(z2,z1,ndum)'; % dummy variable for integration
else
    zdum = zz;
end

switch band_SR
    
    case 0     % test: shortwave radiation are absorbed at the surface
    Iz = zeros(ndum,ntm);
    
    case 2     % 2-band exponential decay of shortwave radiation
    Iz = get_SRz(NSW,zdum,band_SR,waterType);
    
    case 9     % 9-band exponential decay of shortwave radiation
    Iz = get_SRz(NSW,zdum,band_SR,waterType);
end

% depth dependent flux due to the penetrative shortwave radiation
w_theta_r = -(NSW' - Iz)/cp/rho0;
w_b_r     = g*(alpha .* w_theta_r);

%% MOST parameters

% buoyancy forcing
FS.Bf = -(w_b_r + w_b_0);

% surface layer scales
FS.Ustar =  SKF.Ustar';
FS.Sstar = -w_s_0 ./ FS.Ustar / kappa;
FS.Tstar = -(w_theta_r + w_theta_0) ./ FS.Ustar / kappa;
FS.Bstar =  FS.Bf ./ FS.Ustar / kappa;

% Monin-Obukhov length [m]
MOL = (FS.Ustar.^3) ./ (kappa*FS.Bf);

% stability parameter
if numel(zdum) == numel(MOL)
    zeta = zdum' ./ MOL;
else
    zeta = repmat(abs(zdum),1,ntm) ./ MOL;
end

% Append non-turbuelnt fluxes
FS.w_theta_r = w_theta_r;
FS.w_b_r     = w_b_r;

end