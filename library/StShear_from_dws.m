function [duStdz,dvStdz] = StShear_from_dws(WQ,fctr,bw,ziSt,wdir)
%
% StShear_from_dws
%==========================================================================
%
% USAGE:
%  [duStdz,dvStdz] = StShear_from_dws(WQ,fctr,bw,ziSt,wdir)
%
% DESCRIPTION:
%  Compute shear of the Stokes drift from directional wave spectrum.
%  The first band centered at 0.02 Hz is a noise band, hence not used in
%   calculation.
%
% INPUT:
%
%  WQ   - timetable for spectral wave data
%  fctr - bin center frequency
%  bw   - bandwidth for frequency bins
%  ziSt - vertical coordinates for Stokes drift shear
%  wdir - direction of surface wind, clockwise from N [degree, from]
%
% OUTPUT:
%
%  duStdz - Stokes drift shear in x direction
%  dvStdz - Stokes drift shear in y direction
%
% AUTHOR:
%  January 23 2020, Zhihua Zheng                          [ zhihua@uw.edu ]
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
mult1 = 128*pi^5/(g^2);
mult2 = 16*pi^3/g;
hfcoe =  2*pi^3/g;
decoe =  8*pi^2/g;
fc    = fctr(end) + bw(end)/2; % right edge cutoff frequency [Hz]

nz  = length(ziSt);
ntm = length(datm);
nb  = length(bw);

%% Reshaping

df   = repmat(bw,  1,ntm,nz);
f    = repmat(fctr,1,ntm,nz);
z    = repmat(ziSt, 1,ntm,nb);
z    = permute(z,[3 2 1]);

Spec = repmat(spec,1,1,nz);
A1   = repmat(a1,  1,1,nz);
B1   = repmat(b1,  1,1,nz);

%% Shear from resolved Stokes spectrum

decay = exp(decoe*z .* (f.^2));

duStdz_re = -squeeze(mult1*sum(f.^5 .* Spec .* B1 .* decay .* df))';
dvStdz_re = -squeeze(mult1*sum(f.^5 .* Spec .* A1 .* decay .* df))';

%% Shear from analytical high frequency tail (Breivik et al. 2014)

% Note the direction of the waves in the high frequency tail of the 
% spectrum is that of the wind (see Vincent et al. 2019).

% R1 at cutoff frequency
R1c = sqrt(a1(end,:).^2 + b1(end,:).^2);
a1_tail = R1c .* cosd(wdir);
b1_tail = R1c .* sind(wdir);

mu    = -decoe * ziSt;
tailC =  mult2 * fc^5 * ( sqrt(-hfcoe./ziSt) .* erfc(fc*sqrt(mu)) );

duStdz_tail = -tailC * (spec(end,:) .* b1_tail);
dvStdz_tail = -tailC * (spec(end,:) .* a1_tail);

%% Add up

duStdz = duStdz_re + duStdz_tail;
dvStdz = dvStdz_re + dvStdz_tail;

end

