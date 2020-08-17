function read_ocsp_TS(cGotm_in)
%
% Read hourly subsurface TS measurements from the OCSP mooring
% Average the time series in one-hour intervals
%
% cGotm_in: option to create input files for GOTM simulation
%  1 - true
%  0 - false
%
% Averaged data are saved as 'ocsp_prof_1hrPMEL.mat'
%
% Zhihua Zheng, University of Washington, July 15 2019
% =========================================================================

%% General setting

data_dir = './data/';
TS_dir   = [data_dir,'Papa/TS/'];

Tname = fullfile(TS_dir,'t50n145w_hr.cdf');
Sname = fullfile(TS_dir,'s50n145w_hr.cdf');

%% Read variables

time    = ncread(Tname,'time');
lon     = ncread(Tname,'lon');
lat     = ncread(Tname,'lat');

depth_t = ncread(Tname,'depth');
tprof   = ncread(Tname,'T_20');
tprofQ  = ncread(Tname,'QT_5020'); % quality codes

depth_s  = ncread(Sname,'depth');
sprof    = ncread(Sname,'S_41');
sprofQ   = ncread(Sname,'QS_5041'); % quality codes

%% Double precisions

time    = double(time);
depth_t = double(depth_t);
depth_s = double(depth_s);
lon     = double(lon);
lat     = double(lat);

tprof  = double(squeeze(tprof));
sprof  = double(squeeze(sprof));
tprofQ = double(squeeze(tprofQ));
sprofQ = double(squeeze(sprofQ));

%% Pre-processing

% get the reference time
t_unit  = ncreadatt(Tname,'time','units'); % units for 'time'
ref_str = extractAfter(t_unit,'since '); % truncate to get the time string
t_ref = datenum(ref_str,'yyyy-mm-dd HH:MM:SS'); 

time = t_ref + time/24;
dahr = time*24; % in hours
datm = datetime(time,'ConvertFrom','datenum');

% select high quality data
% quality codes ahve category 0, 1, 2, 3, 4, 5
tprof(tprofQ==0 | tprofQ==3 | tprofQ==4 | tprofQ==5) = NaN;
sprof(sprofQ==0 | sprofQ==3 | sprofQ==4 | sprofQ==5) = NaN;

% absolute salinity
SAprof = gsw_SA_from_SP(sprof,depth_s,lon,lat)';

% fill salinity gaps with linear interpolation
sproff   = vert_fill(sprof,depth_s,time);
SAproff  = gsw_SA_from_SP(sproff,depth_s,lon,lat);
SAprofft = interp1(depth_s,SAproff,depth_t,'linear','extrap')';

% interpolated salinity is used to derive potential temperature
PTprof = gsw_pt0_from_t(SAprofft,tprof',depth_t);

%% More variables

g = 9.81;
rho0 = 1025;

% conservative temperature, potential density & buoyancy
CTprof = gsw_CT_from_t(SAprofft,tprof',depth_t);
PDprof = gsw_sigma0(SAprofft,CTprof);
Bprof = -g*(PDprof - rho0)/rho0;

% buoyancy frequency squared [1/s^2]
NSQprof = center_diff(Bprof,-depth_t',2,'mid');

%% Timestable

PROF = timetable(datm,dahr,PTprof,SAprof,PDprof,Bprof,NSQprof);

save([data_dir,'Papa/ocsp_prof_1hrPMEL.mat'],'depth_t','depth_s','PROF');

%% GOTM input files

if cGotm_in
    
gotmdata_root = '~/Documents/GitHub/GOTM/gotmwork/data/';
basecase      = [gotmdata_root,'OCSPapa_20070608-20190616/'];

SPname   = [basecase,'sprof_file.dat'];
TPname   = [basecase,'tprof_file.dat'];
GRIDname = [basecase,'ocsp_cgrid.dat'];

tproff   = vert_fill(tprof,depth_t,time);
PTproff  = gsw_pt0_from_t(SAprofft,tproff',depth_t);
SAprof3  = reshape(SAproff, length(depth_s),1,length(datm));
PTprof3  = reshape(PTproff',length(depth_t),1,length(datm));
datm_str = string(datestr(datm,'yyyy-mm-dd HH:MM:SS'));

write_gotm_ini(SPname,SAprof3,datm_str,-depth_s)
write_gotm_ini(TPname,PTprof3,datm_str,-depth_t)

% z_gotm  = [0:0.5:11 12:1:25 26:2:38 40:3:70 75:5:145]';
z_gotm  = [0:0.25:15 16:1:25 26:2:38 40:3:70 75:5:145]';
zi_gotm = [0; (z_gotm(1:end-1) + z_gotm(2:end))/2; 150];
nLayer  = length(z_gotm);
Ndgotm  = [nLayer; diff(zi_gotm)];

fID = fopen(GRIDname,'w');
fprintf(fID,'%5.3f\n',Ndgotm);
fclose(fID);
end

end