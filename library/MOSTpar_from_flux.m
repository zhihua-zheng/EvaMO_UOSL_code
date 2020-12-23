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
%  idatm    - 1-D column vector, MATLAB datetime for the period inquired
%  zz       - Vertical coordinates for the choosen 2 levels [-, m]
%  casename - Name of the case for which parameters are inquired
%         
% OUTPUT:
%
%  MOL  - Monin-Obukhov length for evaluated depth [m]
%  zeta - Monin-Obukhov stability parameter
%  FS   - Struct contains turbulence scale for velocity (Ustar), 
%         temperature (Tstar), salinity (Sstar), buoyancy (Bstar),
%          buoyancy forcing Bf
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

w_b_0 = SKF.w_b_0';
w_s_0 = SKF.w_s_0';
w_t_0 = SKF.w_t_0';
alpha = SKF.alpha';
NSW   = SKF.nsw;

%% Constants

ntm = length(w_b_0);

g     = 9.81;
rho0  = 1025;
cp    = 3985;
kappa = 0.4;

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
    
    case 0 % test: shortwave radiation are absorbed at the surface
    Iz = zeros(ndum,ntm);
    
    case 2 % 2-band exponential decay of shortwave radiation
    Iz = get_SRz(NSW,zdum,band_SR,waterType);
    
    case 9 % 9-band exponential decay of shortwave radiation
    Iz = get_SRz(NSW,zdum,band_SR,waterType);
end

% depth dependent flux due to the penetrative shortwave radiation
w_t_r = -(NSW' - Iz)/cp/rho0;
w_b_r = g*(alpha .* w_t_r);

%% MOST parameters

% buoyancy forcing
FS.Bf = -(w_b_r + w_b_0);

% surface layer scales
FS.Ustar =  SKF.Ustar';
FS.Sstar = -w_s_0 ./ FS.Ustar / kappa;
FS.Tstar = -(w_t_r + w_t_0) ./ FS.Ustar / kappa;
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
FS.w_t_r = w_t_r;
FS.w_b_r = w_b_r;

end