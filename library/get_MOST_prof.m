function PTmost = get_MOST_prof(z,zdum,sst,Tstar,zeta,fop)
%
% get_MOST_prof
%==========================================================================
%
% USAGE:
%  PTmost = get_MOST_prof(z,sst,Tstar,zeta,fop)
%
% DESCRIPTION:
%  Compute integrated profile from surface to depth based on surface layer
%    scaling
%
% INPUT:
%
%  z     - 1-D column vector (z), vertical coordinates for output profiles,
%          top down [-, m]
%  zdum  - 1-D column vector (z), dummy vertical coordinates for 
%          integration, [-, m]
%  sst   - 1-D row vector(t), surface values for each profile
%  Tstar - 2-D matrix (z,t), surface layer temperature scale
%  zeta  - 2-D matrix (z,t), Monin-Obukhov stability parameter, |z|/MOL
%  fop   - formulation option for surface layer scaling law
%          'Wall'   - neutral law of wall
%          'LMD94'  - Large et al. 1994
%          'Kansas' - Businger et al. 1971
%
% OUTPUT:
%
%  PTmost - 2-D matrix (z,t), integrated profiles
%
% AUTHOR:
%  November 27 2019, Zhihua Zheng                          [ zhihua@uw.edu ]
%==========================================================================

[~,ntm] = size(zeta);
phi_dum = get_emp_phi(zeta,fop);
intgrd  = Tstar ./ abs(zdum) .* phi_dum;

cdPTmost = cumtrapz(zdum,intgrd);

% label the integrated profile as NaN when forcing is not availabe
for j = 1:ntm    
    if sum(~isnan(cdPTmost(2:end,j))) == 0
        cdPTmost(1,j) = NaN;
    end
end

% add surface values
aPTmost = cdPTmost + sst;

% Pick the integrated values at depths of interest
Iout = ismember(zdum,z);
PTmost = aPTmost(Iout,:);

end