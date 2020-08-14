function IntPhiX = get_MOST_Delta(Xstar,MOL,zz,fop)
%
% get_MOST_Delta
%==========================================================================
%
% USAGE:
%  IntPhiX = get_MOST_Delta(Xstar,MOL,zz,fop)
%
% DESCRIPTION:
%  Use the integral approach to estimate the difference of mean quantity
%  at different levels based on Monin-Obukhov scaling with formulation fop
%
% INPUT:
%
%  Xstar - the tubulent scale for the quantity X
%  MOL   - Monin-Obukhov length [m]
%  zz    - vertical coordinates for the choosen 2 levels [-, m]
%  fop   - Monin-Obukhov universal function formulation option
%          'LMD94'  for Large et al. 1994
%          'Kansas' for Businger et al. 1971
%         
% OUTPUT:
%
%  IntPhiX - difference of quantity X at two different levels
%
% AUTHOR:
%  July 24 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Pre-setting

ndum = size(Xstar,1);

z1   = zz(1);
z2   = zz(2);
zdum = linspace(z2,z1,101)'; % dummy variable for integration
zeta = abs(zdum) ./ MOL;

%% Compute empirical similarity function (scalar) & integration

[phiXdum,~] = get_emp_phi(zeta,fop);

if ndum == 1
    Xstar = repmat(Xstar,101,1);
end

intgrd  = phiXdum .* (Xstar ./ abs(zdum));
IntPhiX = trapz(zdum,intgrd);

end