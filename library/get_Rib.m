function Ribprof = get_Rib(Bprof,Ustar,z)
%
% get_Rib
%==========================================================================
%
% USAGE:
%  Ribprof = get_Rib(Bprof,Ustar,z)
%
% DESCRIPTION:
%  Compute nondimensional bulk Richardson number for different levels
%  Rib = N^2 / S^2 = [(Bs - Bh)/(zs - zh)]/[2*Ustar/kappa/(zs + zh)]^2
%
% INPUT:
%
%  Bprof - [z,t] 2-D matrix, time series of buoyancy profile [m^2/s]
%  Ustar - [1,t] 1-D vector, waterside friction velocity [m/s]
%  z     - [z,1] 1-D vector, vertical coordinates [-, m]
%     
% OUTPUT:
%
%  Ribprof - bulk Richardson number at different levels
%
% AUTHOR:
%  August 8 2019, Zhihua Zheng                            [ zhihua@uw.edu ]
%==========================================================================

%% Presetting

kappa = 0.4;
nz  = length(z);
ntm = length(Ustar);

%% The law of the wall approximation for velocity shear squared [1/s^2]

del_z = z(1) - z(2:end);
S  = repmat(Ustar,nz-1,1) ./ repmat(kappa*del_z,1,ntm) .* ...
     repmat(log(z(2:end)/z(1)),1,ntm); % consider roughness length?
S2 = S.^2;

%% Buoyancy frequency squared [1/s^2]

del_B = repmat(Bprof(1,:),nz-1,1) - Bprof(2:end,:);
N2 = del_B ./ repmat(del_z,1,ntm);

%% Bulk Richardson number

Ribprof = N2 ./ S2;

end