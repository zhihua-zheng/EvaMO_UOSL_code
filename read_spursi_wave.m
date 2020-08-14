%% read_spursi_wave
%
% Read directional wave spectrum data from the SPURS-I mooring
% Average the time series in one-hour and three-hours intervals
%
% Averaged data are saved as 'spursi_wave_1hr.mat' and
%  'spursi_wave_3hr.mat', respectively
%
% Zhihua Zheng, University of Washington, July 30 2019
%
% Updates:
%  - Compute vetical shear of the Stokes drift from the directional wave 
%    spectrum, January 23 2020
% =========================================================================

%% General setting

root_dir = '~/GDrive/UW/Research/Data/';
Wave_dir = [root_dir,'SPURSI/Wave/'];
Met_dir  = [root_dir,'SPURSI/Met/'];

Mname    = fullfile(Met_dir, 'SPURS_2012_D_M_1hr.nc');
PROFname = fullfile(root_dir,'SPURSI/spursi_prof_1hrBox.mat');

%% Load data

% frequency bins for the wave spectrum
% 1st bin is for noise, last 3 are dummy 0's
load([Wave_dir,'freq'],'fctr','wave_bw');

DWSfiles = dir([Wave_dir,'/41061_*.mat']);
n_DWS    = length(DWSfiles);
n_f      = length(wave_bw);

wave_time = nan(n_DWS,1);
wave_spec = nan(n_DWS,n_f); % spectral energy density (m^2/hertz)
wave_a1   = nan(n_DWS,n_f);
wave_b1   = nan(n_DWS,n_f);
wave_mdir = nan(n_DWS,n_f);

for i = 1:n_DWS %3761
    
    DWS = load([Wave_dir,DWSfiles(i).name]);
    
    wave_time(i) = DWS.mday; % matlab datenumber
    
    % wave mean direction [degree, from], clockwise from N
    wave_mdir(i,:) = DWS.alpha1;
    
    % nondirectional spectrum
    wave_spec(i,:) = DWS.c11;
    
    % convert to Longuet-Higgins coefficients
    wave_a1(i,:)   = DWS.R1 .* cosd(DWS.alpha1);
    wave_b1(i,:)   = DWS.R1 .* sind(DWS.alpha1);
end

% remove dummy bins
fctr    = fctr(1:end-3);
wave_bw = wave_bw(1:end-3);
wave_spec = wave_spec(:,1:end-3);
wave_a1   = wave_a1(:,1:end-3);
wave_b1   = wave_b1(:,1:end-3);
wave_mdir = wave_mdir(:,1:end-3);

dahr = wave_time*24; % in hours
datm = datetime(wave_time,'ConvertFrom','datenum');

%% Compute the Stokes drift velocity and its vertical shear

wq = timetable(datm,wave_spec,wave_a1,wave_b1);

% construct vertical coordinates
load(PROFname,'depth_t');
zSt     = [flip(-depth_t(1:28)); (-0.6:0.2:0)'];
zSt_mid = (zSt(1:end-1) + zSt(2:end))/2;
ziSt    = union(zSt(1:end-1),zSt_mid);

% load wind measurements
uwind = ncread(Mname,'UWND');
vwind = ncread(Mname,'VWND');
wdir  = atan2d(vwind,uwind); % counterclockwise from E [degree, to]
wdir  = 180+(90-wdir);

% interpolate the wind direction to the time of wave measurements
wind_time = ncread(Mname,'TIME');
wind_time = wind_time + datenum(1950,1,1,0,0,0);
wdir = interp1_ang(wind_time,wdir,wave_time)';

[uSt,vSt] = Stokes_from_dws(wq,fctr,wave_bw,zSt,wdir);
[duSt_dz,dvSt_dz] = StShear_from_dws(wq,fctr,wave_bw,ziSt,wdir);

% wave displacement variance from integral of the nondirectional spectrum
m0 = sum(wave_spec(:,2:end)' .* wave_bw(2:end)); % [m^2]
Hs = 4*sqrt(m0); % significant wave height [m]

%% Timetable

Hs  = Hs';
uSt = uSt';
vSt = vSt';
duSt_dz = duSt_dz';
dvSt_dz = dvSt_dz';

sd = timetable(datm,dahr,uSt,vSt,duSt_dz,dvSt_dz,Hs);

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

save([root_dir,'SPURSI/spursi_wave_1hr.mat'],'SD','zSt');

%% 3-Hour average

% bin edges and center
tlower = dateshift(datm(1),  'start','day');
tupper = dateshift(datm(end),'end',  'day');
datm3E = (tlower:hours(3):tupper)';
datm_3hr = datm3E + hours(1) + minutes(30);

% use 'mean' as there are some gaps in time series (e.g. Jan-31-0000 2013)
SD3 = retime(sd,datm3E,'mean');
SD3.datm = datm_3hr;
SD3.dahr = datenum(datm_3hr)*24; % in hours

% eliminate bins with more/less than 3 data points for 3-hour average
SD = SD3(5:3060,:);

save([root_dir,'SPURSI/spursi_wave_3hr.mat'],'SD','zSt');
