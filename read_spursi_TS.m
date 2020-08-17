function read_spursi_TS
%
% Read subsurface TS measurements from the SPURS-I mooring
% Average the time series in one-hour intervals
%
% Averaged data are saved as 'spursi_prof_1hrBox.mat'
%
% Zhihua Zheng, culy 17 2019
% =========================================================================

%% General setting

data_dir = './data/';
Met_dir  = [data_dir,'SPURSI/Met/'];
TS_dir   = [data_dir,'SPURSI/TS/'];

TSname = fullfile(TS_dir,'SPURS_2012_D_TS.nc'); % 5-min TS profiles

% add surface TS measurements from the file for meteorological data
Mname = fullfile(Met_dir,'SPURS_2012_D_M.nc'); % 1-min resolution

%% Read variables 

% TIME in NC files: days since 1950-01-01T00:00:00Z

time = ncread(TSname,'TIME');
lat  = ncread(TSname,'LATITUDE');
lon  = ncread(TSname,'LONGITUDE');

depth_t = ncread(TSname,'DEPTH') + 0.115; % correction upon recovery
tprof   = ncread(TSname,'TEMP');
depth_s = ncread(TSname,'DEPTH') + 0.115; % correction upon recovery
sprof   = ncread(TSname,'PSAL');
% instM   = ncread(TSname,'INST_MODEL')';

sss = ncread(Mname,'PSAL');
sst = ncread(Mname,'TEMP');
tss = ncread(Mname,'TIME');

% depth of surface TS measurements
depth_ss = ncread(Mname,'DEPTH');

%% Adjust timestamp

time = datenum(1950,1,1,0,0,0) + time;
dahr = time*24; % in hours
datm = datetime(time,'ConvertFrom','datenum'); % dt = 5 min

tss  = datenum(1950,1,1,0,0,0) + tss;
dss  = datetime(tss,'ConvertFrom','datenum'); % dt = 1 min

%% Pre-processing

surfTS = timetable(dss,sst,sss);
datmE  = datm - minutes(2); % bins for 5-min average
SurfTS = retime(surfTS,datmE,@mean);
SurfTS.dss = datm;

% eliminate bins with less than 5 data points for 5-min average
SurfTS(1,:)      = {NaN NaN};
SurfTS(109605,:) = {NaN NaN};

depth_s = [depth_ss; depth_s];
depth_t = [depth_ss; depth_t];
sprof   = [SurfTS.sss'; sprof]; % append sss
tprof   = [SurfTS.sst'; tprof]; % append sst

% truncate time series
Islc  = datm <= datetime(2013,9,1,5,25,0);
sprof = sprof(:,Islc);
tprof = tprof(:,Islc);
time  = time(Islc);
dahr  = dahr(Islc);
datm  = datm(Islc);

% fill salinity gaps with linear interpolation & matrix transposition
sproff = vert_fill(sprof,depth_s,time)';
tprof  = tprof';
sprof  = sprof';

% absolute salinity
SAproff = gsw_SA_from_SP(sproff,depth_s,lon,lat);
SAprof  = gsw_SA_from_SP(sprof, depth_s,lon,lat);

% interpolated salinity is used to derive potential temperature
PTprof  = gsw_pt0_from_t(SAproff,tprof,depth_t);

% aggregate all profiles
Prof = timetable(datm,dahr,PTprof,SAprof,SAproff,tprof,sprof);

%% Hourly average

% bin edges and center
tlower = dateshift(datm(1),  'start','hour');
tupper = dateshift(datm(end),'end',  'hour');
datmE  = (tlower:hours(1):tupper)';
datm_hr = datmE + minutes(30);

PROF1 = retime(Prof,datmE,@mean);
PROF1.datm = datm_hr;
PROF1.dahr = datenum(datm_hr)*24; % in hours

% eliminate bins with more/less than 12 data points for hourly average
PROF = PROF1(2:end-2,:);

%% More variables in hourly timetable

g = 9.81;
rho0 = 1025;

% conservative temperature, potential density & buoyancy
CTprof = gsw_CT_from_t(PROF.SAproff,PROF.tprof,depth_t);
PDprof = gsw_sigma0(PROF.SAproff,CTprof);
Bprof  = -g*(PDprof - rho0)/rho0;

% buoyancy frequency squared [1/s^2]
NSQprof = center_diff(Bprof,-depth_t,2,'mid');

% append potential density, buoyancy and N2
PROF = addvars(PROF,PDprof,Bprof,NSQprof);

save([data_dir,'SPURSI/spursi_prof_1hrBox.mat'],...
     'depth_t','depth_s','PROF');

%% 3-hourly average

% % bin edges and center
% tlower = dateshift(datm(1),  'start','day');
% tupper = dateshift(datm(end),'end',  'day');
% datm3E = (tlower:hours(3):tupper)';
% datm_3hr = datm3E + hours(1) + minutes(30);
% 
% PROF3 = retime(Prof,datm3E,@mean);
% PROF3.datm = datm_3hr;
% PROF3.dahr = datenum(datm_3hr)*24; % in hours
% 
% % eliminate bins with more/less than 36 data points for 3-hour average
% PROF = PROF3(8:2817,:);

%% More variables in 3-hourly timetable

% CTprof = gsw_CT_from_t(PROF.SAproff,PROF.tprof,depth_t);
% PDprof = gsw_sigma0(PROF.SAproff,CTprof);
% Bprof  = -g*(PDprof - rho0)/rho0;
% 
% % buoyancy frequency squared [1/s^2]
% NSQprof = center_diff(Bprof,-depth_t,2,'mid');
% 
% % append potential density, buoyancy and N2
% PROF = addvars(PROF,PDprof,Bprof,NSQprof);
% 
% save([root_dir,'SPURSI/spursi_prof_3hrBox.mat'],...
%      'depth_t','depth_s','PROF');

end