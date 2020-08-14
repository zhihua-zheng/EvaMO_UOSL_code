function [mld,Ilast] = get_mld(A,z,flag,avar)
%
% get_mld
%==========================================================================
%
% USAGE:
%  mld = get_mld(A,z,flag,z_ref)
%
% DESCRIPTION:
%  Compute the mixed layer depth from the density, temperature, or bulk 
%   Richardson number profile.
%
% INPUT:
%
%  A - 2-D matrix (z,t) containing density, temperature or Rib profiles
%  z - 1-D column vector, vertical coordinates for profiles, bottom up
%      [-, m]
%  flag - 1 density criteria (de Boyer Montégut et al. 2004)
%         2 temperature criteria (de Boyer Montégut et al. 2004)
%         3 bulk Richardson number criteria
%  avar - reference depth (-1/-10) or, critical bulk Richardson number
%     
% OUTPUT:
%
%  mld - mixed layer depth [+, m]
%
% AUTHOR:
%  October 29 2018, Zhihua Zheng                          [ zhihua@uw.edu ]
%==========================================================================

[nz,ntm] = size(A);
        
switch flag
    case 1       
        if A(1) > 100
            sigma = A - 1000;
        else
            sigma = A;
        end
        
        % Get reference values
        sigma_ref = zeros(1,ntm);
        for j = 1:ntm
            goodi = ~isnan(sigma(:,j));
            if sum(goodi) < 2
                sigma_ref(j) = NaN;
            else
                sigma_ref(j) = interp1(z(goodi),sigma(goodi,j),avar,'linear','extrap');
            end
        end

        % MLD is where the density is 0.03 kg/m^3 higher than
        % a surface reference value at 10 m (see de Boyer Montégut et al.
        % 2004), or at 1 m (considering shallow surface layer)   
        dval = sigma_ref + 0.03;
        where_mld = sigma - repmat(dval,nz,1);
    
    case 2
        temp = A;

        % Get reference values
        temp_ref = zeros(1,ntm);
        for j = 1:ntm
            goodi = ~isnan(temp(:,j));
            if sum(goodi) < 2
                temp_ref(j) = NaN;
            else
                temp_ref(j) = interp1(z(goodi),temp(goodi,j),avar,'linear','extrap');
            end
        end

        % MLD as the shallowest depth where the temp. is 0.2 C colder than
        % a surface reference value at 10 m (de Boyer Montégut et al. 2004),
        % or at 1 m (considering shallow surface layer)   
        dval = temp_ref - 0.2;
        where_mld = repmat(dval,nz,1) - temp;
        
    case 3
        Rib = A;
        Ric = avar;
        where_mld = Rib - Ric;
end
  
%% interpolate to find where the zero is

mld   = zeros(ntm,1);
Ilast = zeros(ntm,1);

for j = 1:ntm

    mldj = where_mld(:,j);
    tmp  = mldj(~isnan(mldj));

    if length(tmp) < 2 % not enough points to interpolate
        mld(j) = NaN;
    else
        tmpz = z(~isnan(mldj));
        mldD = find(tmp>0,1,'last');
        if isempty(mldD) % mixed all the way to bottom
            Ilast(j) = 1;
            mld(j) = -z(1);
            continue
        else
            Ilast(j) = mldD+1; % index of last level in mixed layer
        end
        
        if mldD == length(tmp)
            mld(j) = NaN; % mld is shallower than upmost level, unresolvable here
        else
            mld(j) = -interp1(tmp(mldD:mldD+1),tmpz(mldD:mldD+1),0);
        end
    end
end

end
