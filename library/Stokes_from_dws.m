function [uSt,vSt] = Stokes_from_dws(WQ,fctr,bw,zSt)
%
% Stokes_from_dws
%==========================================================================
%
% USAGE:
%  [uSt,vSt] = Stokes_from_dws(WQ,fctr,bw,zSt,wdir)
%
% DESCRIPTION:
%  Compute the Stokes drift velocity from directional wave spectrum.
%
% INPUT:
%
%  WQ   - timetable for spectral wave data
%  fctr - bin center frequency
%  bw   - bandwidth for frequency bins
%  zSt  - [z,1] vertical coordinates for the Stokes drift
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
a1   = WQ.wave_a1';
b1   = WQ.wave_b1';

%% Constants

g     = 9.81;
multi = 16*pi^3/g;
decoe =  8*pi^2/g;
fc    = fctr(end) + bw(end)/2; % right edge cutoff frequency [Hz]

nz  = length(zSt);
ntm = length(datm);
nb  = length(bw);

%% Reshaping

df = repmat(bw,  1,ntm,nz);
f  = repmat(fctr,1,ntm,nz);
z  = repmat(zSt, 1,ntm,nb);
z  = permute(z,[3 2 1]);

A1 = repmat(a1,1,1,nz);
B1 = repmat(b1,1,1,nz);

%% Resolved Stokes spectrum

decay = exp(decoe*z .* (f.^2));

uSt_re = -squeeze(multi*sum(f.^3 * pi .* B1 .* decay .* df))';
vSt_re = -squeeze(multi*sum(f.^3 * pi .* A1 .* decay .* df))';

%% Append analytical high frequency tail (Breivik et al. 2014)

a1_tail = a1(end,:);
b1_tail = b1(end,:);

mu    = -decoe * zSt;
tailC =  multi * fc^4 * ( exp(-mu*fc^2) - ...
         fc*sqrt(mu*pi) .* erfc(fc*sqrt(mu)) );

uSt_tail = -tailC * (pi * b1_tail);
vSt_tail = -tailC * (pi * a1_tail);

%% Add up

uSt = uSt_re + uSt_tail;
vSt = vSt_re + vSt_tail;

end