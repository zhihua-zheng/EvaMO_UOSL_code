function [phiS,phiM,TKE] = get_SMCse_phi_nol(zeta,mov)
%
% get_SMCse_phi_nol
%==========================================================================
%
% USAGE:
%  [phiS,phiM,TKE] = get_SMCse_phi_nol(zeta,mov)
%
% DESCRIPTION:
%  Compute Monin-Obukhov similarity functions for scalars and momentum
%   from the super-equilibrium (SE) limit (MY level 2, local balance of 
%   turbulence production and dissipation) of Second Moment Closure (SMC).
%  Length scale is assumed to be kappa|z|
%   
% INPUT:
%
%  zeta - Monin-Obukhov stability parameter, |z|/MOL
%  mov  - Model version
%        'MY1982' for Mellor-Yamada 1982 model 
%        'KC1994' for Kansa-Clayson 1994 model
%     
% OUTPUT:
%
%  phiS - Similarity function for scalars  from the SE limit of SMC
%  phiM - Similarity function for momentum from the SE limit of SMC
%
% AUTHOR:
%  December 13 2019, Zhihua Zheng                         [ zhihua@uw.edu ]
%==========================================================================

%% Constants

if strcmp(mov,'MY1982')

A1 = 0.92;
A2 = 0.74;
B1 = 16.6;
B2 = 10.1;
C1 = 0.08;
C2 = 0;
C3 = 0;

elseif strcmp(mov,'KC1994')
    
A1 = 0.92;
A2 = 0.74;
B1 = 16.6;
B2 = 10.1;
C1 = 0.08;
C2 = 0.7;
C3 = 0.2;

elseif strcmp(mov,'KC2004')

A1 = 0.92;
A2 = 0.74;
B1 = 16.6;
B2 = 10.1;
C1 = 0.08;
C2 = 0.7;
C3 = 0.2;

else
    disp('model version not available!');
    return
end

%% Combined constants

r  = 1/3 - 2*A1/B1;
R1 = r - C1;

Z3 = (6*A1 + B2*(1-C3));
Z2 =  3*A2;
Z1 = (6*A1 + 3*A2*(1-C2));

S2 = 1/3/A2;
S1 = 1/3/A1;

%% Dimensionaless functions

nzet = numel(zeta);
phiM = nan(size(zeta));
phiS = nan(size(zeta));

fun = @SMCse_phiEq_nol;
options = optimoptions('fsolve',...'Algorithm','Levenberg-Marquardt',...
                       'Display','off');

for i = 1:nzet
    
    zt = zeta(i);
    
    % change starting point to avoid singularity
    if zt > 0
        p0(1:2) = zt*[5, 5] + 1;
    else
        p0(1:2) = [1, 1];
    end
    p0(3) = (B1*(p0(1)-zt))^(1/3);
    
    P = fsolve(fun,p0,options);
    if P(3)<0
        disp('unrealistic solution occured!');
        return
    end
    
    phiM(i) = P(1);
    phiS(i) = P(2);
    
    % normalized TKE components
    TKE.uu(i) = P(3)^2*r + 6*A1*P(1)/P(3);
    TKE.vv(i) = P(3)^2*r;
    TKE.ww(i) = P(3)^2*r - 6*A1*zt/P(3);
end

%% Construct system of equations

function F = SMCse_phiEq_nol(p)

%==========================================================================
% p(1) : \phi_m
% p(2) : \phi_h
% p(3) : q^*
%==========================================================================

F(1) = p(1)*( R1 - Z1*zt/(p(3)^3) ) - p(2)*Z2*zt/(p(3)^3) - S1/p(3);
   
F(2) = p(2)*( r - Z3*zt/(p(3)^3) ) - S2/p(3);

F(3) = B1*( p(1) - zt ) - p(3)^3;

end

end