%% read_spursi_flux
%
% Read computed air-sea fluxes from the SPURS-I mooring
% Average the time series in one-hour intervals
%
% Averaged data are saved as 'spursi_flux_1hrUOP.mat'
%
% Zhihua Zheng, University of Washington, July 17 2019
% =========================================================================

%% General setting

data_dir = './data/';
Met_dir  = [data_dir,'SPURSI/Met/'];

% hourly meterological variables and fluxes
Mname  = fullfile(Met_dir,'SPURS_2012_D_M_1hr.nc');
SFname = fullfile(Met_dir,'SPURS_2012_D_F_1hr.nc');

%% Read variables

% TIME in NC files: days since 1950-01-01T00:00:00Z

time = ncread(SFname,'TIME');
lon  = ncread(Mname,'LONGITUDE');
lat  = ncread(Mname,'LATITUDE');

tau  = ncread(SFname,'TAUMAG'); % [N/m^2]
taud = ncread(SFname,'TAUDIR'); % [degree, to], clockwise from N
nsw  = ncread(SFname,'QS');     % [W/m^2] into the ocean
nlw  = ncread(SFname,'QL');     % [W/m^2] into the ocean
hlb  = ncread(SFname,'QH');     % [W/m^2] into the ocean
hsb  = ncread(SFname,'QB');     % [W/m^2] into the ocean
tsk  = ncread(SFname,'TSKIN');  % [C] skin temperature

uwin = ncread(Mname,'UWND');    % [m/s] eastward  wind speed
vwin = ncread(Mname,'VWND');    % [m/s] northward wind speed
zwin = ncread(Mname,'HEIGHT_WND'); % [m] height of wind measurements
airt = ncread(Mname,'AIRT');    % [C]   air temperature
sst  = ncread(Mname,'TEMP');    % [C]   surface temperature
sss  = ncread(Mname,'PSAL');    % [PSU] surface salinity
ssp  = ncread(Mname,'DEPTH');   % [m]   depth of surface TS measurements
rain = ncread(Mname,'RAIN');    % [mm/hr]

% absolute salinity
ssSA = gsw_SA_from_SP(sss,ssp,lon,lat); % [g/kg]

% heat of evaporation for seawater
Le = gsw_latentheat_evap_t(ssSA,tsk); % [J/kg]

rho_fw = 1000;
evap   = -hlb ./ Le / rho_fw; % [m/s], mostly positive
evap   = evap*1000*3600;      % [mm/hr]

tau_x  = tau .* sind(taud);
tau_y  = tau .* cosd(taud);
wspd   = sqrt(uwin.^2 + vwin.^2);
u10    = spshfttc(wspd,zwin,10,airt);

%% Adjust timestamp

time = datenum(1950,1,1,0,0,0) + time;
dahr = time*24; % in hours
datm = datetime(time,'ConvertFrom','datenum');
sf = timetable(datm,dahr,tau_x,tau_y,tau,u10,...
               nsw,nlw,hlb,hsb,evap,rain,sst,sss);

%% Hourly average

% bin edges and center
tlower = dateshift(datm(1),  'start','hour');
tupper = dateshift(datm(end),'end',  'hour');
datmE  = (tlower:hours(1):tupper)';
datm_hr = datmE + minutes(30);

SF1 = retime(sf,datmE,@mean);
SF1.datm = datm_hr;
SF1.dahr = datenum(datm_hr)*24; % in hours

% eliminate empty bin
SF = SF1(1:end-1,:);

save([data_dir,'SPURSI/spursi_flux_1hrUOP.mat'],'SF');
clear

%% 3-hourly average

% % bin edges and center
% tlower = dateshift(datm(1),  'start','day');
% tupper = dateshift(datm(end),'end',  'day');
% datm3E = (tlower:hours(3):tupper)';
% datm_3hr = datm3E + hours(1) + minutes(30);
% 
% SF3 = retime(sf,datm3E,@mean);
% SF3.datm = datm_3hr;
% SF3.dahr = datenum(datm_3hr)*24; % in hours
% 
% % eliminate bins with more/less than 3 data points for 3-hour average
% SF = SF3(8:3051,:);
% 
% save([data_dir,'SPURSI/spursi_flux_3hrUOP.mat'],'SF');
