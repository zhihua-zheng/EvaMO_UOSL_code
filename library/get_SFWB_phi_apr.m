function [phiS,phiM,TKE] = get_SFWB_phi_apr(zeta,xi,alpha_b,mop)
%
% get_SFWB_phi_apr
%==========================================================================
%
% USAGE:
%  [phiS,phiM,TKE] = get_SFWB_phi_apr(zeta,xi,alpha_b,mop)
%
% DESCRIPTION:
%  Compute Monin-Obukhov similarity functions for scalars and momentum
%   from Second Moment Closure including surface TKE flux from surface 
%   wave breaking (Craig and Banner 1994 + buoyancy flux). 
%  q^3 is set by the approximate solution in equation (18).
%   
% INPUT:
%
%  zeta - Monin-Obukhov stability parameter, |z|/MOL
%  xi   - Normalized vertical distance, |z|/z_0
%  mop  - Model option:
%         1 model SFWB as turbulent transport  of TKE
%         2 model SFWB as turbulent transports of TKE and TKE-components
%         3 model SFWB as pressure  transport  of TKE (Janssen 2012)
%
% OUTPUT:
%
%  phiS - Similarity function for scalars
%  phiM - Similarity function for momentum
%  TKE  - Three components of the TKE, normalized by u_*^2
%
% AUTHOR:
%  April 8 2020, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Parsing inputs

nZet = numel(zeta);
nXi  = numel(xi);
nAlB = numel(alpha_b);

if nAlB == 1
    alpha_b = repmat(alpha_b,size(zeta));
elseif nAlB ~= nZet
    disp('Unable to interpret parameter alpha_b ...')
    disp('Please check input dimension!')
    return
end

if nXi == 1
    xi = repmat(xi,size(zeta));
elseif nZet == 1
    zeta = repmat(zeta,size(xi));    
elseif nXi ~= nZet
    disp('Unable to interpret parameter xi ...')
    disp('Please check input dimension!')
    return
end

%% Single constants

A1 = 0.92;
A2 = 0.74;
B1 = 16.6;
B2 = 10.1;
C1 = 0.08;
C2 = 0.7;
C3 = 0.2;

Sq = 0.2;
kappa = 0.4;

%% Combined constants

c  = alpha_b/2*sqrt(3*B1/Sq);
n  = sqrt(3/B1/Sq/(kappa^2));

r  = 1/3 - 2*A1/B1;
R1 = r - C1;

Z3 = (6*A1 + B2*(1-C3));
Z2 =  3*A2;
Z1 = (6*A1 + 3*A2*(1-C2));

X2 =  6*A1/5/B1;
X1 = 18*A1/5/B1;

S2 = 1/3/A2;
S1 = 1/3/A1;

%% Dimensionaless functions

phiM   = nan(size(zeta));
phiS   = nan(size(zeta));
TKE.uu = nan(size(zeta));
TKE.vv = nan(size(zeta));
TKE.ww = nan(size(zeta));
TKE.qq = nan(size(zeta));

switch mop
    case 1
        fun = @SFWB_phiEq1_nol;    
    case 2
        fun = @SFWB_phiEq2_nol;        
    case 3
        fun = @SFWB_phiEq3_nol;
end

options = optimoptions('fsolve',...'Algorithm','Levenberg-Marquardt',...
                       'Display','off');

for j = 1:numel(zeta) 
    
    zt = zeta(j);
    
    if mop == 3
        eb  = xi(j);
        alB = alpha_b(j)*kappa;
    else
        eb = c(j)*xi(j)^(-n);
    end
    
    aPar = [zt eb];
    
    if sum(isnan(aPar)) == 0
        
        % change starting point to avoid singularity
        if zt > 0
            p0(1:2) = zt*[1, 1] + 5;
        else
            p0(1:2) = [2.5, 2.5];
        end
        p0(3) = (B1*5 + eb)^(1/3);
        
        P = fsolve(fun,p0,options);
        if P(3)<0
            disp('unrealistic solution occured!');
            return
        end

        phiM(j) = P(1);
        phiS(j) = P(2);
        
        switch mop
            case 1
            % normalized TKE components for mop = 1
            TKE.uu(j) = P(3)^2*r + 6*A1/P(3)*P(1);
            TKE.vv(j) = P(3)^2*r;
            TKE.ww(j) = P(3)^2*r - 6*A1/P(3)*zt;
            TKE.qq(j) = P(3)^2;

            case 2
            % normalized TKE components for mop = 2
            TKE.uu(j) = P(3)^2*r + 6*A1/P(3)*P(1) + X2/P(3)*eb;
            TKE.vv(j) = P(3)^2*r + X2/P(3)*eb;
            TKE.ww(j) = P(3)^2*r - 6*A1/P(3)*zt + X1/P(3)*eb;
            TKE.qq(j) = P(3)^2;
        end
        
    end
end

%% Construct system of equations

function F = SFWB_phiEq1_nol(p)

%==========================================================================
% p(1) : \phi_m^x
% p(2) : \phi_h
% p(3) : q^*
%==========================================================================

F(1) = p(1)*( R1 - Z1*zt/(p(3)^3) ) - p(2)*Z2*zt/(p(3)^3) - S1/p(3);
   
F(2) = p(2)*( r - Z3*zt/(p(3)^3) ) - S2/p(3);

F(3) = B1*( p(1) - zt ) + eb - p(3)^3;

end

%% Construct system of equations

% this version includes turbulent diffusion in TKE-component equations

function F = SFWB_phiEq2_nol(p)

%==========================================================================
% p(1) : \phi_m^x
% p(2) : \phi_h
% p(3) : q^*
%==========================================================================

F(1) = p(1)*( R1 - (Z1*zt - X1*eb)/(p(3)^3) ) - ...
       p(2)*Z2*zt/(p(3)^3) - S1/p(3);
   
F(2) = p(2)*( r - (Z3*zt - X1*eb)/(p(3)^3) ) - S2/p(3);

F(3) = B1*( p(1) - zt ) + eb - p(3)^3;

end

%% Construct system of equations

% this version includes pressure diffusion in TKE equation

function F = SFWB_phiEq3_nol(p)

%==========================================================================
% p(1) : \phi_m^x
% p(2) : \phi_h
% p(3) : q^*
%==========================================================================

F(1) = p(1)*( R1 - Z1*zt/(p(3)^3) ) - p(2)*Z2*zt/(p(3)^3) - S1/p(3);
   
F(2) = p(2)*( r - Z3*zt/(p(3)^3) ) - S2/p(3);

F(3) = B1*( p(1)*(1-exp(-4*eb))^2 - zt + alB*eb*exp(-eb) ) - p(3)^3;

end

end