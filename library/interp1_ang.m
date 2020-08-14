function angd_q = interp1_ang(x,angd,xq)
%
% interp1_ang
%==========================================================================
%
% USAGE:
%  angd_q = interp1_ang(x,angd,xq)
%
% DESCRIPTION:
%  Linearly interpolate angles to specific query points
%
% INPUT:
%
%  x    - 1-D vector, sample points of angle 
%  angd - 1-D vector, sampled angle [degree]
%  xq   - 1-D vector, query points of angle
% 
% OUTPUT:
%
%  angd_q - 1-D vector, interpolated angle at xq [degree]
%
% AUTHOR:
%  June 26 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

ang   = deg2rad(angd); % convert to radian
d_ang = angdiff(ang);
d_x   = diff(x);

nq = numel(xq);
ang_q = nan(size(xq));

for j = 1:nq
    
    Il = find(x-xq(j)<0,1,'last');
    
    dq_x   = xq(j) - x(Il);
    dq_ang = dq_x/d_x(Il-1)*d_ang(Il-1);
    ang_q(j) = ang(Il) + dq_ang;
end

angd_q = rad2deg(ang_q); % convert to degree

end