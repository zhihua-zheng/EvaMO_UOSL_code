function [uSt,vSt] = Stokes_from_dws(WQ,fctr,bw,zSt,wdir)
%
% Stokes_from_dws
%==========================================================================
%
% USAGE:
%  [uSt,vSt] = Stokes_from_dws(WQ,fctr,bw,zSt,wdir)
%
% DESCRIPTION:
%  Compute the Stokes drift velocity from directional wave spectrum.
%  The first band centered at 0.02 Hz is a noise band, hence not used in
%   calculation.
%
% INPUT:
%
%  WQ   - timetable for spectral wave data
%  fctr - bin center frequency
%  bw   - bandwidth for frequency bins
%  zSt  - [z,1] vertical coordinates for the Stokes drift
%  wdir - [1,t] direction of surface wind, clockwise from N [degree, from]
%
% OUTPUT:
%
%  uSt - Stokes drift profile in x direction
%  vSt - Stokes drift profile in y direction
%
% AUTHOR:
%  July 31 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%% Loading

datm = WQ.datm;
spec = WQ.wave_spec(:,2:end)';
a1   = WQ.wave_a1(:,2:end)';
b1   = WQ.wave_b1(:,2:end)';

fctr = fctr(2:end);
bw   = bw(2:end);

%% Constants

g     = 9.81;
multi = 16*pi^3/g;
decoe =  8*pi^2/g;
fc    = fctr(end) + bw(end)/2; % right edge cutoff frequency [Hz]

nz  = length(zSt);
ntm = length(datm);
nb  = length(bw);

%% Reshaping

df   = repmat(bw,  1,ntm,nz);
f    = repmat(fctr,1,ntm,nz);
z    = repmat(zSt, 1,ntm,nb);
z    = permute(z,[3 2 1]);

Spec = repmat(spec,1,1,nz);
A1   = repmat(a1,  1,1,nz);
B1   = repmat(b1,  1,1,nz);

%% Resolved Stokes spectrum

decay = exp(decoe*z .* (f.^2));

uSt_re = -squeeze(multi*sum(f.^3 .* Spec .* B1 .* decay .* df))';
vSt_re = -squeeze(multi*sum(f.^3 .* Spec .* A1 .* decay .* df))';

%% Append analytical high frequency tail (Breivik et al. 2014)

% Note the direction of the waves in the high frequency tail of the 
% spectrum is that of the wind. Vincent et al. 2019

% R1 at cutoff frequency
R1c = sqrt(a1(end,:).^2 + b1(end,:).^2);
a1_tail = R1c .* cosd(wdir);
b1_tail = R1c .* sind(wdir);

mu    = -decoe * zSt;
tailC =  multi * fc^4 * ( exp(-mu*fc^2) - ...
         fc*sqrt(mu*pi) .* erfc(fc*sqrt(mu)) );

uSt_tail = -tailC * (spec(end,:) .* b1_tail);
vSt_tail = -tailC * (spec(end,:) .* a1_tail);

%% Add up

uSt = uSt_re + uSt_tail;
vSt = vSt_re + vSt_tail;

end

