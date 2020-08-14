function read_ocsp_wave(cGotm_in)
%
% Read directional wave spectrum data from the OCSP mooring
% Average the time series in one-hour and three-hours intervals
%
% cGotm_in: option to create input files for GOTM simulation
%  1 - true
%  0 - false
%
% Averaged data are saved as 'ocsp_wave_1hr.mat' and
%  'ocsp_wave_3hr.mat', respectively
%
% Zhihua Zheng, University of Washington, August 13 2018
%
% Updates:
%  - Compute vetical shear of the Stokes drift from the directional wave 
%    spectrum, Feburary 4 2020
% =========================================================================

%% General setting

root_dir = '~/GDrive/UW/Research/Data/';
Wave_dir = [root_dir,'Papa/Wave'];
Met_dir  = [root_dir,'Papa/Met'];

DWSname  = fullfile(Wave_dir,'166p1_historic.nc');
WINDname = fullfile(Met_dir, 'w50n145w_hr.cdf');
PROFname = fullfile(root_dir,'Papa/ocsp_prof_1hrPMEL.mat');

%% Read variables 

wave_time = double(ncread(DWSname,'waveTime'));

% band center frequency [Hz]
fctr = double(ncread(DWSname,'waveFrequency'));

% frequency bandwidth [Hz]
wave_bw = ncread(DWSname,'waveBandwidth');

% nondirectional spectral energy density [m^2/Hz]
wave_spec = (ncread(DWSname,'waveEnergyDensity'))'; % equivalent to C_11

% Longuet-Higgins coefficients
wave_a1 = (ncread(DWSname,'waveA1Value'))';
wave_b1 = (ncread(DWSname,'waveB1Value'))';

% significant wave height [m]
Hs = ncread(DWSname,'waveHs');

%% Adjust timestamp

% get the reference time
t_unit  = ncreadatt(DWSname,'waveTime','units'); % unit for 'time'
ref_str = extractAfter(t_unit,'since '); % truncate to get the time string
t_ref   = datenum(ref_str,'yyyy-mm-dd HH:MM:SS'); 

wave_time = t_ref + double(wave_time)/3600/24;
dahr = wave_time*24; % in hours
datm = datetime(wave_time,'ConvertFrom','datenum');

%% the Stokes drift velocity and its vertical shear

wq = timetable(datm,wave_spec,wave_a1,wave_b1);

% construct vertical coordinates
load(PROFname,'depth_t');
dmid    = (depth_t(1:end-1) + depth_t(2:end))/2;
depth   = union(dmid,depth_t);
zSt     = [flip(-depth(1:31)); (-0.8:0.2:0)'];
zSt_mid = (zSt(1:end-1) + zSt(2:end))/2;
ziSt    = union(zSt(1:end-1),zSt_mid);

nzSt  = length(zSt);
nziSt = length(ziSt);
ntm   = length(datm);

uSt = nan(nzSt,ntm);
vSt = nan(nzSt,ntm);
duSt_dz = nan(nziSt,ntm);
dvSt_dz = nan(nziSt,ntm);

% load wind direction [degree, to], clockwise from N 
wind_time = ncread(WINDname,'time');
wind_time = double(wind_time)/24 + datenum(2007,6,8,4,0,0);
wind_dir  = double(squeeze(ncread(WINDname,'WD_410'))); 
Qwind_dir = double(squeeze(ncread(WINDname,'QWD_5410')));
wind_dir(Qwind_dir==0) = NaN;
wind_dir = wind_dir + 180; % clockwise from N [degree, from]

% interpolate the wind direction to the time of wave measurements
wdir = interp1_ang(wind_time,wind_dir,wave_time)';

% split matrix into chunks due to the MATLAB size limit
nChunk = fix(ntm/10000);
jt = (0:nChunk)'*10000;
jt(1) = 1;
jt(end) = ntm + 1;

for j = 1:nChunk
    
    jl = jt(j);
    jr = jt(j+1) - 1;
    [uSt(:,jl:jr),...
     vSt(:,jl:jr)] = Stokes_from_dws(wq(jl:jr,:),fctr,...
                       wave_bw,zSt,wdir(jl:jr));
    [duSt_dz(:,jl:jr),...
     dvSt_dz(:,jl:jr)] = StShear_from_dws(wq(jl:jr,:),fctr,...
                           wave_bw,ziSt,wdir(jl:jr));
end

%% Timetable

uSt = uSt';
vSt = vSt';
duSt_dz = duSt_dz';
dvSt_dz = dvSt_dz';

sd = timetable(datm,dahr,uSt,vSt,duSt_dz,dvSt_dz,Hs);

%% Hourly average

% bin edges and center
tlower  = dateshift(datm(1),  'start','hour');
tupper  = dateshift(datm(end),'end',  'hour');
datm_hr = (tlower:hours(1):tupper)';
datmE   = datm_hr - minutes(30);

SD1 = retime(sd,datmE,'mean'); % hourly average
SD1.datm = datm_hr;
SD1.dahr = datenum(datm_hr)*24; % in hours

% remove empty bins
SD = SD1(2:end-1,:);

save([root_dir,'Papa/ocsp_wave_1hr.mat'],'SD','zSt');

%% 3-Hour average

% bin edges and center
tlower   = dateshift(datm(1),  'start','day');
tupper   = dateshift(datm(end),'end',  'day');
datm_3hr = (tlower:hours(3):tupper)';
datm3E   = datm_3hr - hours(1) - minutes(30);

SD3 = retime(sd,datm3E,'mean');
SD3.datm = datm_3hr;
SD3.dahr = datenum(datm_3hr)*24; % in hours

% eliminate bins with more/less than 3 data points for 3-hour average
SD = SD3(7:end-8,:);

save([root_dir,'Papa/ocsp_wave_3hr.mat'],'SD','zSt');

%% GOTM input file

if cGotm_in
    
gotmdata_root = '~/Documents/GitHub/GOTM/gotmwork/data/';
basecase      = [gotmdata_root,'OCSPapa_20070608-20190616/'];
WSPname       = [basecase,'spec_file.dat']; % wave spectrum file

% band mean direction [degree, from], clockwise from N
wave_mdir = (ncread(DWSname,'waveMeanDirection'))';
% equivalent to theta_1 = atan2(wave_b1,wave_a1)*180/pi + 360; (verified!)

% convert wave direction to counter-clockwise from E, [degree, to]
wave_mdir = 90 - (180 + wave_mdir);

wave_r1 = sqrt(wave_a1.^2 + wave_b1.^2);
xcmp    = cosd(wave_mdir)'; 
ycmp    = sind(wave_mdir)'; 
df      = repmat(wave_bw,1,length(datm));
spec_h2 = wave_spec' .* wave_r1' .* df;

waveInput  = cat(3,spec_h2,xcmp,ycmp);
waveInput3 = permute(waveInput,[1 3 2]);
datm_str   = string(datestr(datm,'yyyy-mm-dd HH:MM:SS'));

write_gotm_ini(WSPname,waveInput3,datm_str,fctr)
gzip(WSPname)
delete(WSPname)
end

end