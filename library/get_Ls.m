function Ls = get_Ls(ziSt,BLD,dUStDw_dz,kappa_l)
%
% get_Ls
%==========================================================================
%
% USAGE:
%  Ls = get_Ls(ziSt,BLD,dUStDw_dz,kappa_l)
%
% DESCRIPTION:
%  Compute the CL vortex forcing length scale in in-homogenenous
%  pressure-strain rate closure (Harcourt 2015)
%
% INPUT:
%
%  ziSt      - 1-D vector (z), vertical coordinates for Stokes shear
%  BLD       - 1-D vector (t), mixed layer depth
%  dUStDw_dz - 2-D matrix (t,z), downwind Stokes shear
%  kappa_l   - assumed ratio of the length scale to |z|
%
% OUTPUT:
%
%  Ls - 1-D vector (1,t), the CL vortex forcing length scale
%
% AUTHOR:
%  January 21 2020, Zhihua Zheng                          [ zhihua@uw.edu ]
%==========================================================================

%% Read relevant variables

StDw_shear = dUStDw_dz';
[nzS,ntm]  = size(StDw_shear);

%% CL vortex forcing length scale Ls

Ls = nan(1,ntm);

% parabolic shaped length scale in the OSBL
l_parab = kappa_l*abs(ziSt).*(ziSt+BLD')./BLD';
l = l_parab;
l(l<0) = NaN;

% shape function for turbulent momentum flux
Gwu = 1 + ziSt ./ BLD';
Gwu(Gwu<0) = 0; % nullify the momentum flux outside boundary layer

% positive CL production weighted length scale
for j = 1:ntm
    jStDw_shear = StDw_shear(:,j);
    jGwu = Gwu(:,j);
    
    if sum(isnan(jStDw_shear)) == 0
        where_PCL = jGwu .* jStDw_shear;
        
        % lelvel below positive CL production
        ibPCL = find(where_PCL<0,1,'last');
        if isempty(ibPCL) % all positive case
            ibPCL = find(where_PCL<1e-4,1,'last');
        end
        iPos = ibPCL + 1;

        if iPos < nzS
            Intd_num = jGwu(iPos:nzS).*jStDw_shear(iPos:nzS).*l(iPos:nzS,j);
            Intd_den = jGwu(iPos:nzS).*jStDw_shear(iPos:nzS);
            
            % note the integration doesn't end exactly at surface
            Num = trapz(ziSt(iPos:nzS),Intd_num);
            Den = trapz(ziSt(iPos:nzS),Intd_den);
            Ls(j) = Num/Den;
        end
    end
end

end