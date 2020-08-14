function A_filled = vert_fill(A,z,t)
%
% vert_fill
%==========================================================================
%
% USAGE:
%  A_filled = vert_fill(A,z,t)
%
% DESCRIPTION:
%  Function to vertically fill the missing data, using linear interpolation
%  as default method.
%
% INPUT:
%  
%  A - 2-D matrix, Quantity to be interpolated
%  z - 1-D vector, row coordinates of A [-, m]
%  t - 1-D vector, column coordinates of A [+]
%
% OUTPUT:
%
%  A_filled - matrix A with gaps filled
%
% AUTHOR:
%  March 12 2019. Zhihua Zheng                            [ zhihua@uw.edu ]
%==========================================================================

%% Horizontally fill the surface and bottom

A_surf = A(1,:);
A_bot  = A(end,:);

if sum(isnan(A_surf)) >= 1 && sum(~isnan(A_surf)) > 0
    A(1,:) = interp1(t(~isnan(A_surf)),A_surf(~isnan(A_surf)),t);
end

if sum(isnan(A_bot)) >= 1 && sum(~isnan(A_bot)) > 0
    A(end,:) = interp1(t(~isnan(A_bot)),A_bot(~isnan(A_bot)),t);
end

A_filled = A;

%% Vertically fill the columns

[~, cols] = find(isnan(A));
col_nan = unique(cols);

for j = 1:length(col_nan)
    
    tmp   = A(:,col_nan(j));
    Igood = ~isnan(tmp);
    
    if sum(Igood) <= 2
        continue
    end
    
    A_filled(:,col_nan(j)) = interp1(z(Igood),tmp(Igood),z); %,'linear','extrap'
end
    
end