function [duStdz,dvStdz] = StShear_from_dws(WQ,fctr,bw,ziSt)
%
% StShear_from_dws
%==========================================================================
%
% USAGE:
%  [duStdz,dvStdz] = StShear_from_dws(WQ,fctr,bw,ziSt,wdir)
%
% DESCRIPTION:
%  Compute shear of the Stokes drift from directional wave spectrum.
%
% INPUT:
%
%  WQ   - timetable for spectral wave data
%  fctr - bin center frequency
%  bw   - bandwidth for frequency bins
%  ziSt - vertical coordinates for Stokes drift shear
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
a1   = WQ.wave_a1';
b1   = WQ.wave_b1';

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

df = repmat(bw,  1,ntm,nz);
f  = repmat(fctr,1,ntm,nz);
z  = repmat(ziSt, 1,ntm,nb);
z  = permute(z,[3 2 1]);

A1 = repmat(a1,1,1,nz);
B1 = repmat(b1,1,1,nz);

%% Shear from resolved Stokes spectrum

decay = exp(decoe*z .* (f.^2));

duStdz_re = -squeeze(mult1*sum(f.^5 * pi .* B1 .* decay .* df))';
dvStdz_re = -squeeze(mult1*sum(f.^5 * pi .* A1 .* decay .* df))';

%% Shear from analytical high frequency tail (Breivik et al. 2014)

a1_tail = a1(end,:);
b1_tail = b1(end,:);

mu    = -decoe * ziSt;
tailC =  mult2 * fc^5 * ( sqrt(-hfcoe./ziSt) .* erfc(fc*sqrt(mu)) );

duStdz_tail = -tailC * (pi * b1_tail);
dvStdz_tail = -tailC * (pi * a1_tail);

%% Add up

duStdz = duStdz_re + duStdz_tail;
dvStdz = dvStdz_re + dvStdz_tail;

end