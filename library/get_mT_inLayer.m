function mTil = get_mT_inLayer(PTprof,depth_t,H,fun)
%
% get_mT_inLayer
%==========================================================================
%
% USAGE:
%  mTil = get_mT_inLayer(PTprof,depth_t,H,fun)
%
% DESCRIPTION:
%  Compute the mean/min/max quantity in a layer with thickness H.
%
% INPUT:
%
%  PTprof  - 2-D matrix (z,t) containing the qunatity of interest
%  depth_t - 1-D column vector with vertical coordinates for profiles, 
%            top down [+, m]
%  H       - 1-D column of the layer thickness [+, m]
%  fun     - function handle (e.g., @nanmean, @nanmin, @nanmax, etc.)
%     
% OUTPUT:
%
%  mTml - mean/min/max value in the layer
%
% AUTHOR:
%  November 27 2019, Zhihua Zheng                         [ zhihua@uw.edu ]
%==========================================================================

[~,ntm] = size(PTprof);
mTil = nan(1,ntm);

for j = 1:ntm
    
    h = H(j);
    if isnan(h)
        mTil(j) = NaN;
    else
        where_ld = depth_t - h;
        Ild = find(where_ld < 0,1,'last');
        if ~isempty(Ild)
            mTil(j) = fun(PTprof(1:Ild,j));
        end
    end
end

end