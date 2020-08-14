function [phiS,phiM] = get_emp_phi(zeta,fop)
%
% get_emp_phi
%==========================================================================
%
% USAGE:
%  [phiS,phiM] = get_emp_phi(zeta,fop)
%
% DESCRIPTION:
%  Compute empirical Monin-Obukhov similarity functions for scalars and
%  momentum according to the formulation option 'fop'.
%
% INPUT:
%
%  zeta - Monin-Obukhov stability parameter, |z|/MOL
%  fop - formulation option
%        'Wall' for neutral law of wall
%        'LMD94' for Large et al. 1994
%        'Kansas' for Businger et al. 1971
%        'Hog88' for Kansas curve corrected by Högström 1988
%        'SPURSI' for fitted function from SPURS-I
%
% OUTPUT:
%
%  phiS - Empirical similarity function for scalars
%  phiM - Empirical similarity function for momentum
%
% AUTHOR:
%  July 19 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Constants

if strcmp(fop,'LMD94')

zetaS = -1;
cS    =  98.96;
aS    = -28.86;

zetaM = -0.2;
cM    =  8.38;
aM    =  1.26;

elseif strcmp(fop,'Kansas') || strcmp(fop,'Hog88') || strcmp(fop,'SPURSI')
    
zetaS = -2;
zetaM = -2;

elseif strcmp(fop,'Wall')

phiS = ones(size(zeta));
phiM = ones(size(zeta));
return

else
    
disp('formulation option not available!');
phiS = nan(size(zeta));
phiM = nan(size(zeta));
return

end

%% Empirical dimensionaless functions

phiS = nan(size(zeta));
phiM = nan(size(zeta));

iStable     = find(zeta >= 0);
iWUnstableS = find(zeta < 0 & zeta >= zetaS);
iWUnstableM = find(zeta < 0 & zeta >= zetaM);
iSUnstableS = find(zeta <= zetaS);
iSUnstableM = find(zeta <= zetaM);

% stable side
if strcmp(fop,'LMD94')
    phiS(iStable) = 1.00 + 5.0*zeta(iStable);
    phiM(iStable) = 1.00 + 5.0*zeta(iStable);
    
elseif strcmp(fop,'Kansas')
    phiS(iStable) = 0.74 + 4.7*zeta(iStable);
    phiM(iStable) = 1.00 + 4.7*zeta(iStable);
    
elseif strcmp(fop,'Hog88')
    phiS(iStable) = 1.00 + 7.8*zeta(iStable); % after correction
    phiM(iStable) = 1.00 + 6.0*zeta(iStable); % after correction

elseif strcmp(fop,'SPURSI')
    phiS(iStable) = 0.4 + 5*zeta(iStable);
    phiM(iStable) = 0.4 + 5*zeta(iStable);

else
    disp('formulation option not available!');
end

% unstable side
if strcmp(fop,'LMD94')
    phiS(iWUnstableS) = (1 - 16*zeta(iWUnstableS)).^(-1/2);
    phiM(iWUnstableM) = (1 - 16*zeta(iWUnstableM)).^(-1/4);
    
    phiS(iSUnstableS) = (aS - cS*zeta(iSUnstableS)).^(-1/3);
    phiM(iSUnstableM) = (aM - cM*zeta(iSUnstableM)).^(-1/3);
    
elseif strcmp(fop,'Kansas')
    phiS(iWUnstableS) = 0.74*(1 -  9*zeta(iWUnstableS)).^(-1/2);
    phiM(iWUnstableM) =      (1 - 15*zeta(iWUnstableM)).^(-1/4);

    phiS(iSUnstableS) = 0.74*(1 -  9*zeta(iSUnstableS)).^(-1/2);
    phiM(iSUnstableM) =      (1 - 15*zeta(iSUnstableM)).^(-1/4);

elseif strcmp(fop,'Hog88')
    phiS(iWUnstableS) = (1 - 12.0*zeta(iWUnstableS)).^(-1/2); % after correction
    phiM(iWUnstableM) = (1 - 19.3*zeta(iWUnstableM)).^(-1/4); % after correction

    phiS(iSUnstableS) = (1 - 12.0*zeta(iSUnstableS)).^(-1/2); % after correction
    phiM(iSUnstableM) = (1 - 19.3*zeta(iSUnstableM)).^(-1/4); % after correction

elseif strcmp(fop,'SPURSI')
    phiS(iWUnstableS) = 0.4*(1 - 3.6*zeta(iWUnstableS)).^(-1/2);
    phiM(iWUnstableM) = 0.4*(1 - 3.6*zeta(iWUnstableM)).^(-1/4);

    phiS(iSUnstableS) = 0.4*(1 - 3.6*zeta(iSUnstableS)).^(-1/2);
    phiM(iSUnstableM) = 0.4*(1 - 3.6*zeta(iSUnstableM)).^(-1/4);

else
    disp('formulation option not available!');
end


end

