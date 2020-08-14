function [phiS,phiM,rs,TKE] = get_H15se_phi(zeta,etaX,etaY,fzS,E6,mop)
%
% get_H15se_phi
%==========================================================================
%
% USAGE:
%  [phiS,phiM,rs,TKE] = get_H15se_phi(zeta,etaX,etaY,fzS,E6,mop)
%
% DESCRIPTION:
%  Compute Monin-Obukhov similarity functions for scalars and momentum
%   from the super-equilibrium (SE) limit (MY level 2, local balance of 
%   turbulence production and dissipation) of Second Moment Closure 
%   including Langmuir Turbulence (SMCLT).
%  Model version: Harcourt 2015.
%  Length scale is set by algebraic equation (25) in the paper.
%   
% INPUT:
%
%  zeta - Monin-Obukhov stability parameter, |z|/MOL
%  etaX - Downwind  Langmuir stability parameter, -(du^s/dz)/(u_*/kappa|z|) 
%  etaY - Crosswind Langmuir stability parameter, -(dv^s/dz)/(u_*/kappa|z|)
%  fzS  - Surface proximity function, 1 + tanh(z/4/Ls)
%         Decay length Ls is computed by function "get_Ls"
%  mop  - Model option:
%         1 wave aligned with wind 
%         2 general wave direction
%         3 general wave direction & transport of vertical TKE (<w'w'>)
%
% OUTPUT:
%
%  phiS - Similarity function for scalars  from the SE limit of SMCLT
%  phiM - Similarity function for momentum from the SE limit of SMCLT
%  rs   - Ratio of turbulent length scale to kappa|z|
%  TKE  - Three components of the TKE, normalized by u_*^2
% 
% AUTHOR:
%  December 14 2019, Zhihua Zheng                         [ zhihua@uw.edu ]
%==========================================================================

%% Single constants

A1  = 0.92;
A2  = 0.74;
B1  = 16.6;
B2  = 10.1;
C1  = 0.08;
C2  = 0.7;
C3  = 0.2;

E1  = 1.8;
E2  = 1.0;
E3  = 1.8;
E4  = 1.33;

Sl  = 0.2;

%% Combined constants

r  = 1/3 - 2*A1/B1;
R1 = r - C1;

Z4 = (6*A1 + B2*(1-C3));
Z2 =  3*A2;
Z1 = (6*A1 + 3*A2*(1-C2));

T6 =  9*A1;
T5 =  6*A1; % different than that in "get_H13se_phi"
T4 =  3*A2*(1-C2);
T3 = (6*A1 + 3*A2);
T2 =  3*A1;
T1 = 12*A1;

S2 = 1/3/A2;
S1 = 1/3/A1;

%% Dimensionless functions

nZet = numel(zeta);
nEtX = numel(etaX);
nEtY = numel(etaY);
nFzS = numel(fzS);
    
if nEtX == 1
    etaX = repmat(etaX,size(zeta));
elseif nEtX ~= nZet  
    disp('Unable to interpret Langmuir stability parameter \eta^x ...')
    disp('Please check input dimension!')
    return
end
   
if nEtY == 1
    etaY = repmat(etaY,size(zeta));
elseif nEtY ~= nZet  
    disp('Unable to interpret Langmuir stability parameter \eta^y ...')
    disp('Please check input dimension!')
    return
end

if nFzS == 1
    fzS = repmat(fzS,size(zeta));
    
elseif nFzS ~= nZet  
    disp('Unable to interpret surface proximity function ...')
    disp('Please check input dimension!')
    return
end

phiM = nan(size(zeta));
phiS = nan(size(zeta));
rs   = nan(size(zeta));

switch mop
    case 1
        fun  = @H15se_phiEq1;
        tfac = 1-fzS;
    
    case 2
        fun  = @H15se_phiEq2;
        tfac = 1-fzS;
    
    case 3
        fun  = @H15se_phiEq3;
        tfac = 1-fzS/2;
end

kappa = 0.4;
options = optimoptions('fsolve',...'Algorithm','Levenberg-Marquardt',...
                       'Display','off');

for j = 1:nZet   
    
    zt  = zeta(j);
    etX = etaX(j);
    etY = etaY(j);
    tX  = etX*tfac(j);
    tY  = etY*tfac(j);
    aPar = [zt etX tX tY];
    
    if sum(isnan(aPar)) == 0
    
        p0(4:5) = [tY, 1];

        if mop == 3
            bp = zt+tX;
        else
            bp = zt+etX;
        end
        
        % change starting point to avoid singularity
        if bp > 0
            p0(1:2) = bp*[1, 1] + 5;
        else
            p0(1:2) = [2.5, 2.5];
        end
        p0(3) = (B1*5)^(1/3);

        P = fsolve(fun,p0,options);
        if P(3)<0
            disp('unrealistic solution occured!');
            return
        end

        phiM(j) = P(1);
        phiS(j) = P(2);
        rs(j)   = P(5);
        
        % normalized TKE components (\gamma = pi/2) for mop = 1 or 2
        TKE.uu(j) = P(3)^2*r + T5*P(5)/P(3)*P(1);
        TKE.vv(j) = P(3)^2*r - T5*P(5)/P(3)*(etX - tX);
        TKE.ww(j) = P(3)^2*r - T5*P(5)/P(3)*(zt + tX);
    end
end

%% Construct system of equations for misaligned wind - waves

function F = H15se_phiEq2(p)

%==========================================================================
% p(1) : \phi_m^x
% p(2) : \phi_h
% p(3) : q^*
% p(4) : \phi_m^y
% p(5) : l^*
%==========================================================================

F(1) = p(1)*( R1 - (Z1*zt + T1*tX)*p(5)/(p(3)^3) ) - ...
       p(2)*Z2*zt*p(5)/(p(3)^3) - S1/p(5)/p(3) - ...
       p(4)*T2*tY*p(5)/(p(3)^3) - tX*r;

F(2) = p(2)*( r - (Z4*zt + T3*tX)*p(5)/(p(3)^3) ) - ...
       p(1)*T4*tX*p(5)/(p(3)^3) - S2/p(5)/p(3) - ...
       p(4)*T4*tY*p(5)/(p(3)^3);

F(3) = B1*p(5)*( p(1) - zt - etX ) - p(3)^3;

F(4) = p(4)*( R1 - (Z1*zt + T6*tX)*p(5)/(p(3)^3) ) - tY*r;

F(5) = p(5)*( E1*p(1) - E3*zt - E6*etX ) - ...
       E2*(p(3)^3)*( 1 + E4*p(5)^2 )/B1 + ...
       Sl*(p(3)^3)*(p(5)*kappa)^2;
end

%% Construct system of equations for misaligned wind - waves

% this version includes transport of vertical TKE

function F = H15se_phiEq3(p)

%==========================================================================
% p(1) : \phi_m^x
% p(2) : \phi_h
% p(3) : q^*
% p(4) : \phi_m^y
% p(5) : l^*
%==========================================================================

F(1) = p(1)*( R1 - (Z1*zt + T1*tX + T5*(tX-etX))*p(5)/(p(3)^3) ) - ...
       p(2)*Z2*zt*p(5)/(p(3)^3) - S1/p(5)/p(3) - ...
       p(4)*T2*tY*p(5)/(p(3)^3) - tX*r;

F(2) = p(2)*( r - (Z4*zt + T3*tX + T5*(tX-etX))*p(5)/(p(3)^3) ) - ...
       p(1)*T4*tX*p(5)/(p(3)^3) - S2/p(5)/p(3) - ...
       p(4)*T4*tY*p(5)/(p(3)^3);

F(3) = B1*p(5)*( p(1) - zt - tX ) - p(3)^3;

F(4) = p(4)*( R1 - (Z1*zt + T6*tX + T5*(tX-etX))*p(5)/(p(3)^3) ) - tY*r;

F(5) = p(5)*( E1*p(1) - E3*zt - E6*tX ) - ...
       E2*(p(3)^3)*( 1 + E4*p(5)^2 )/B1 + ...
       Sl*(p(3)^3)*(p(5)*kappa)^2;
end

%% Construct system of equations for aligned wind - waves

function F = H15se_phiEq1(p)

%==========================================================================
% p(1) : \phi_m
% p(2) : \phi_h
% p(3) : q^*
% p(4) : \phi_m^y
% p(5) : l^*
%==========================================================================

F(1) = p(1)*( R1 - (Z1*zt + T1*tX)*p(5)/(p(3)^3) ) - ...
       p(2)*Z2*zt*p(5)/(p(3)^3) - S1/p(5)/p(3) - tX*r;

F(2) = p(2)*( r - (Z4*zt + T3*tX)*p(5)/(p(3)^3) ) - ...
       p(1)*T4*tX*p(5)/(p(3)^3) - S2/p(5)/p(3);

F(3) = B1*p(5)*( p(1) - zt - etX ) - p(3)^3;

F(4) = p(4);

F(5) = p(5)*( E1*p(1) - E3*zt - E6*etX ) - ...
       E2*(p(3)^3)*( 1 + E4*p(5)^2 )/B1 + ...
       Sl*(p(3)^3)*(p(5)*kappa)^2;
end

end