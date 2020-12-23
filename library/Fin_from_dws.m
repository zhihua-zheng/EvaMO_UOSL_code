function Fin = Fin_from_dws(WQ,fctr,bw,ws,wdir,fop)
%
% Fin_from_dws
%==========================================================================
%
% USAGE:
%  Fin = Fin_from_dws(WQ,fctr,bw,ws,wdir,fop)
%
% DESCRIPTION:
%  Compute the rate of wind energy input accroding to Plant 1982, or 
%   Donelan and Pierson 1987.
%
% INPUT:
%
%  WQ   - timetable for spectral wave data
%  fctr - bin center frequency
%  bw   - bandwidth for frequency bins
%  ws   - [1,t] wind stress, or [f,t] wind speed at height pi/k
%  wdir - [1,t] direction of wind, clockwise from N [degree, from]
%  fop  - formula option:
%         'P82':  Plant 1982
%         'DP87': Donelan and Pierson 1987
%
% OUTPUT:
%
%  Fin - the rate of energy input from wind
%
% AUTHOR:
%  Sep 6 2020, Zhihua Zheng                               [ zhihua@uw.edu ]
%==========================================================================

%% Loading

datm = WQ.datm;
a0   = WQ.wave_spec'/pi;
a1   = WQ.wave_a1';
b1   = WQ.wave_b1';
a2   = WQ.wave_a2';
b2   = WQ.wave_b2';

%% Constants

rhow = 1025;
rhoa = 1.225;
g    = 9.81;
fc   = fctr(end) + bw(end)/2; % right edge cutoff frequency [Hz]
ntm  = length(datm);

%% Reshaping

df = repmat(bw,  1,ntm);
f  = repmat(fctr,1,ntm);

%% Resolved spectrum

switch fop
    case 'P82'
        multi = 0.04*8*pi^3/g/rhoa.*ws;
        mu = pi*(cosd(wdir).*a1 + sind(wdir).*b1);
        Fin_re = multi .* sum(f.^3 .* mu .* df);
    
    case 'DP87'
        c = g/2/pi./f;
        multi = 0.194*8*pi^3*rhoa/rhow/g;
        mu = pi*(ws.^2.*(a0 + cosd(2*wdir).*a2 + sind(2*wdir).*b2)/2 + ...
                 a0.*c.^2 - ...
                 2*ws.*c.*(cosd(wdir).*a1 + sind(wdir).*b1));
        Fin_re = multi * sum(f.^3 .* mu .* df);
end

%% Append analytical high frequency tail (Breivik et al. 2014)

a1_tail = a1(end,:);
b1_tail = b1(end,:);

tailC = multi * fc^4;

switch fop
    case 'P82'
        muC = pi*(cosd(wdir).*a1_tail + sind(wdir).*b1_tail);
        Fin_tail = tailC .* muC;
        
    case 'DP87'
        muC = pi*(ws(end,:).^2.*(a0(end,:) + cosd(2*wdir).*a2(end,:) + ...
                                 sind(2*wdir).*b2(end,:))/2 + ...
                  a0(end,:)*c(end)^2 - ...
                  2*ws(end,:)*c(end).*(cosd(wdir).*a1_tail + ...
                                       sind(wdir).*b1_tail));
        Fin_tail = tailC * muC;
end

%% Add up

Fin = Fin_re + Fin_tail;

end