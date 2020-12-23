%% read_spursi_wave
%
% Read directional wave spectrum data from the SPURS-I mooring
% Average the time series in one-hour intervals
%
% Averaged data are saved as 'spursi_wave_1hr.mat'
%
% Zhihua Zheng, University of Washington, July 30 2019
%
% Updates:
%  - Compute vetical shear of the Stokes drift from the directional wave 
%    spectrum, January 23 2020
% =========================================================================

%% General setting

data_dir = './data/';
Wave_dir = [data_dir,'SPURSI/Wave/'];
Met_dir  = [data_dir,'SPURSI/Met/'];

Mname    = fullfile(Met_dir, 'SPURS_2012_D_M_1hr.nc');
Fname    = fullfile(Met_dir, 'SPURS_2012_D_F_1hr.nc');
PROFname = fullfile(data_dir,'SPURSI/spursi_prof_1hrBox.mat');

%% Load data

% frequency bins for the wave spectrum
% 1st bin is for noise, last 3 are dummy 0's
load([Wave_dir,'freq'],'fctr','wave_bw');

DWSfiles = dir([Wave_dir,'/41061_*.mat']);
n_DWS    = length(DWSfiles);
n_f      = length(wave_bw);

wave_time = nan(n_DWS,1);
wave_spec = nan(n_DWS,n_f); % spectral energy density (m^2/hertz)
wave_fp   = nan(n_DWS,1);
wave_a1   = nan(n_DWS,n_f);
wave_b1   = nan(n_DWS,n_f);
wave_a2   = nan(n_DWS,n_f);
wave_b2   = nan(n_DWS,n_f);
wave_md   = nan(n_DWS,n_f);
wave_pd   = nan(n_DWS,n_f);

for i = 1:n_DWS %3761
    
    DWS = load([Wave_dir,DWSfiles(i).name]);
    
    wave_time(i) = DWS.mday; % matlab datenumber
    
    % wave mean and principal direction [degree, from], clockwise from N
    wave_md(i,:) = DWS.alpha1;
    wave_pd(i,:) = DWS.alpha2;
    
    % nondirectional spectrum
    wave_spec(i,:) = DWS.c11;
    
    % peak frequency
    [~,imax] = max(DWS.c11);
    wave_fp(i) = fctr(imax);
    
    % convert to Longuet-Higgins coefficients
    % wave_a0(i,:) = DWS.c11/pi;
    wave_a1(i,:) = DWS.R1 .* (DWS.c11/pi) .* cosd(DWS.alpha1);
    wave_b1(i,:) = DWS.R1 .* (DWS.c11/pi) .* sind(DWS.alpha1);
    wave_a2(i,:) = DWS.R2 .* (DWS.c11/pi) .* cosd(DWS.alpha2);
    wave_b2(i,:) = DWS.R2 .* (DWS.c11/pi) .* sind(DWS.alpha2);    
end

% remove dummy bins
fctr = fctr(2:end-3);
wave_bw = wave_bw(2:end-3);
wave_a1 = wave_a1(:,2:end-3);
wave_b1 = wave_b1(:,2:end-3);
wave_a2 = wave_a2(:,2:end-3);
wave_b2 = wave_b2(:,2:end-3);
wave_md = wave_md(:,2:end-3);
wave_pd = wave_pd(:,2:end-3);
wave_spec = wave_spec(:,2:end-3);

dahr = wave_time*24; % in hours
datm = datetime(wave_time,'ConvertFrom','datenum');

% wave_a1n = wave_a1./(wave_spec/pi);
% wave_b1n = wave_b1./(wave_spec/pi);
% wave_r1  = sqrt(wave_b1n.^2 + wave_a1n.^2);
% 
% figure('position',[125 139 904 607]);
% for i = 1:16
% 
% subplot(4,4,i)
% is = 24*(i-1)+1;
% ie = 24*i;
% 
% plot(fctr,wave_r1(is:ie,:),'color',[.8 .8 .8])
% ylim([0 1]); hold on; grid on
% plot(fctr,mean(wave_r1(is:ie,:)),'linewidth',2.5,'color',rgb('coral'))
% 
% end
% [~,hsupY1] = suplabel('$R_1$','y');
% [~,hsupY2] = suplabel('frequency','x');
% set(hsupY1,'fontsize',22,'Interpreter','latex','Units','Normalized')
% set(hsupY2,'fontsize',22,'Interpreter','latex','Units','Normalized')
% set(hsupY1,'Position',hsupY1.Position + [0.04 0 0]);
% set(hsupY2,'Position',hsupY2.Position + [0 0.02 0]);

%% Stokes drift velocity and its vertical shear and wind energy input

wq = timetable(datm,wave_spec,wave_a1,wave_b1,wave_a2,wave_b2);

% construct vertical coordinates
load(PROFname,'depth_t');
zSt     = [flip(-depth_t(1:28)); (-0.6:0.2:0)'];
zSt_mid = (zSt(1:end-1) + zSt(2:end))/2;
ziSt    = union(zSt(1:end-1),zSt_mid);

% load meteorological measurements
wstr     = ncread(Fname,'TAUMAG');
wstr_dir = ncread(Fname,'TAUDIR');
wstr_dir = wstr_dir-180;
airt  = ncread(Mname,'AIRT');
uwind = ncread(Mname,'UWND');
vwind = ncread(Mname,'VWND');
zwind = ncread(Mname,'HEIGHT_WND');
wspd  = sqrt(uwind.^2 + vwind.^2);
wdir  = atan2d(vwind,uwind); % counterclockwise from E [degree, to]
wdir  = 180+(90-wdir);

% interpolate meteorolofical measurements to the time of wave measurements
wind_time = ncread(Mname,'TIME');
wind_time = wind_time + datenum(1950,1,1,0,0,0);
wdir = interp1_ang(wind_time,wdir,wave_time)';
wspd = interp1(wind_time,wspd,wave_time)';
airt = interp1(wind_time,airt,wave_time)';

tau_time = ncread(Fname,'TIME');
tau_time = tau_time + datenum(1950,1,1,0,0,0);
wstr = interp1(tau_time,wstr,wave_time)';

% adjust wind speed to the height of half wave length
g = 9.81;
zHWL = g/4/pi./(fctr.^2);
wspdHWL = spshfttc(wspd,zwind,zHWL,airt);

F82 = Fin_from_dws(wq,fctr,wave_bw,wstr,   wdir,'P82');
F87 = Fin_from_dws(wq,fctr,wave_bw,wspdHWL,wdir,'DP87');
[uSt,vSt] = Stokes_from_dws(wq,fctr,wave_bw,zSt);
[duSt_dz,dvSt_dz] = StShear_from_dws(wq,fctr,wave_bw,ziSt);

% wave displacement variance from the integral of nondirectional spectra
m0 = sum(wave_spec(:,2:end)' .* wave_bw(2:end)); % [m^2]
Hs = 4*sqrt(m0); % significant wave height [m]

% phase speed at the peak [m/s]
Cp = g/2/pi./wave_fp;

%% Timetable

F82 = F82';
F87 = F87';
Hs  = Hs';
uSt = uSt';
vSt = vSt';
duSt_dz = duSt_dz';
dvSt_dz = dvSt_dz';

sd = timetable(datm,dahr,uSt,vSt,duSt_dz,dvSt_dz,Hs,Cp,F82,F87);

%% Hourly average

% bin edges and center
tlower = dateshift(datm(1),  'start','hour');
tupper = dateshift(datm(end),'end',  'hour');
datmE  = (tlower:hours(1):tupper)';
datm_hr = datmE + minutes(30);

% use 'mean' as there are some gaps in the middle (e.g., Jan-31-0000 2013)
SD1 = retime(sd,datmE,'mean');
SD1.datm = datm_hr;
SD1.dahr = datenum(datm_hr)*24; % in hours

% remove empty bin
SD = SD1(1:end-1,:);

save([data_dir,'SPURSI/spursi_wave_1hr.mat'],'SD','zSt');
clear

%% 3-Hour average

% % bin edges and center
% tlower = dateshift(datm(1),  'start','day');
% tupper = dateshift(datm(end),'end',  'day');
% datm3E = (tlower:hours(3):tupper)';
% datm_3hr = datm3E + hours(1) + minutes(30);
% 
% % use 'mean' as there are some gaps in time series (e.g. Jan-31-0000 2013)
% SD3 = retime(sd,datm3E,'mean');
% SD3.datm = datm_3hr;
% SD3.dahr = datenum(datm_3hr)*24; % in hours
% 
% % eliminate bins with more/less than 3 data points for 3-hour average
% SD = SD3(5:3060,:);
% 
% save([data_dir,'SPURSI/spursi_wave_3hr.mat'],'SD','zSt');
