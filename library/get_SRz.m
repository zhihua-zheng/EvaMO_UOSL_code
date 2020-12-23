function Iz = get_SRz(I0,z,bn,Jtype)
%
% get_SRz
%==========================================================================
%
% USAGE:
%  Iz = get_SRz(I0,z,bn,Jtype)
%
% DESCRIPTION:
%  Compute solar radiation at depth z, according to Paulson and Simpson 
%  (1977) parameterization, with modification suggested by Soloviev and 
%  Schlüssel (1966).
%
% INPUT:
%
%  I0 - 1D column vector of surface solar irradiance [W/m^2]
%  z - 1D vector of the depth at which the solar radiation is to be 
%    evaluated [-,m]
%  bn - band number of the model (2 or 9)
%  Jtype - values 1 through 5 corresponds to Jerlov's classification of
%    water types I, IA, IB, II, and III, respectively
% 
% OUTPUT:
%
%  Iz - magnitude of solar radiation at level z
%
% AUTHOR:
%  June 23 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Constants

% spectral weighting coefficients
r9 = [0.2370 0.3600 0.1790 0.0870 0.0800 0.0246 0.0250 0.0070 0.0004];

% irrdadiance absorption coefficients [1/m]
mu9 = ones(size(r9));
mu9(2:end) = [4.405e-1 3.175e1 1.825e2 1.201e3 7.937e3 3.195e3 ...
              1.279e4  6.944e4];

switch Jtype
    case 1
        mu9(1) = 0.066;
        r2     = [0.58 0.42];
        mu2    = 1 ./ [0.35 23];
    case 2
        mu9(1) = 0.076;
        r2     = [0.62 0.38];
        mu2    = 1 ./ [0.60 20];
    case 3
        mu9(1) = 0.088;
        r2     = [0.67 0.33];
        mu2    = 1 ./ [1.0 17];
    case 4
        mu9(1) = 0.132;
        r2     = [0.77 0.23];
        mu2    = 1 ./ [1.5 14];
    case 5
        mu9(1) = 0.382;
        r2     = [0.78 0.22];
        mu2    = 1 ./ [1.4 7.9];
end 

%% Band number 

nz = length(z);

if bn == 9
    
    R  = repmat(r9, nz,1);
    MU = repmat(mu9,nz,1);
    
elseif bn == 2
    
    R  = repmat(r2, nz,1);
    MU = repmat(mu2,nz,1);
    
else
    print('Band number not supported yet!')
end

Z  = repmat(z(:),1,bn);
fr = sum(R .* exp(MU .* Z),2);

% Note possible bug when zdum has same length as time vector

if length(fr) == length(I0) % one depth for one time point
    Iz = fr' .* I0'; % 1-D row vector of Iz at depth |z|
else % multiple depths for one time point
    Iz = fr * I0'; % 2-D (z,t) matrix of Iz at each depth |z|
end
    
end

