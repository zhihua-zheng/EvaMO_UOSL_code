function read_ocsp_wave(cGotm_in)
%
% Read directional wave spectrum data from the OCSP mooring
% Average the time series in one-hour intervals
%
% cGotm_in: option to create input files for GOTM simulation
%  1 - true
%  0 - false
%
% Averaged data are saved as 'ocsp_wave_1hr.mat'
%
% Zhihua Zheng, University of Washington, August 13 2018
%
% Updates:
%  - Compute vetical shear of the Stokes drift from the directional wave 
%    spectrum, Feburary 4 2020
% =========================================================================

%% General setting

data_dir = './data/';
Wave_dir = [data_dir,'Papa/Wave'];
Met_dir  = [data_dir,'Papa/Met'];

DWSname  = fullfile(Wave_dir,'166p1_historic.nc');
WINDname = fullfile(Met_dir, 'w50n145w_hr.cdf');
TAUname  = fullfile(Met_dir, 'tau50n145w_hr.cdf');
ARTname  = fullfile(Met_dir, 'airt50n145w_hr.cdf');
PROFname = fullfile(data_dir,'Papa/ocsp_prof_1hrPMEL.mat');

%% Read variables 

wave_time = double(ncread(DWSname,'waveTime'));

% band center frequency [Hz]
fctr = double(ncread(DWSname,'waveFrequency'));

% frequency bandwidth [Hz]
wave_bw = ncread(DWSname,'waveBandwidth');

% nondirectional spectral energy density [m^2/Hz]
wave_spec = (ncread(DWSname,'waveEnergyDensity'))'; % equivalent to C_11
wave_a0   = wave_spec/pi;

% normalized Longuet-Higgins coefficients
wave_a1n = (ncread(DWSname,'waveA1Value'))';
wave_b1n = (ncread(DWSname,'waveB1Value'))';
wave_a2n = (ncread(DWSname,'waveA2Value'))';
wave_b2n = (ncread(DWSname,'waveB2Value'))';

wave_a1 = wave_a1n.* wave_a0;
wave_b1 = wave_b1n.* wave_a0;
wave_a2 = wave_a2n.* wave_a0;
wave_b2 = wave_b2n.* wave_a0;

% wave_pd = atan2(wave_b2,wave_a2)/2;
% wave_md = (ncread(DWSname,'waveMeanDirection'))';
% wave_sp = (ncread(DWSname,'waveSpread'))';

% significant wave height [m]
Hs = ncread(DWSname,'waveHs');

% peak wave period [s]
Tp = ncread(DWSname,'waveTp');

% phase speed at the peak [m/s]
g  = 9.81;
Cp = g/2/pi.*Tp;

% wave_r1 = sqrt(wave_b1n.^2 + wave_a1n.^2);
% figure('position',[125 139 904 607]);
% for i = 1:16
% 
% subplot(4,4,i)
% is = 48*(i-1)+1;
% ie = 48*i;
% 
% plot(fctr,wave_r1(is:ie,:),'color',[.8 .8 .8])
% ylim([0 1]); hold on; grid on
% plot(fctr,mean(wave_r1(is:ie,:)),'linewidth',2.5,'color',rgb('azure'))
% 
% end
% [~,hsupY1] = suplabel('$R_1$','y');
% [~,hsupY2] = suplabel('frequency','x');
% set(hsupY1,'fontsize',22,'Interpreter','latex','Units','Normalized')
% set(hsupY2,'fontsize',22,'Interpreter','latex','Units','Normalized')
% set(hsupY1,'Position',hsupY1.Position + [0.04 0 0]);
% set(hsupY2,'Position',hsupY2.Position + [0 0.02 0]);

%% Adjust timestamp

% get the reference time
t_unit  = ncreadatt(DWSname,'waveTime','units'); % unit for 'time'
ref_str = extractAfter(t_unit,'since '); % truncate to get the time string
t_ref   = datenum(ref_str,'yyyy-mm-dd HH:MM:SS'); 

wave_time = t_ref + double(wave_time)/3600/24;
dahr = wave_time*24; % in hours
datm = datetime(wave_time,'ConvertFrom','datenum');

%% Stokes drift velocity and its vertical shear and wind energy input

wq = timetable(datm,wave_spec,wave_a1,wave_b1,wave_a2,wave_b2);

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

% load meteorological measurements
tau_time  = ncread(TAUname,'time');
tau_time  = double(tau_time)/24 + datenum(2007,6,8,5,0,0);
wind_str  = double(squeeze(ncread(TAUname,'TAU_440')));
Qwind_str = double(squeeze(ncread(TAUname,'QTAU_5440')));

airt_time = ncread(ARTname,'time');
airt_time = double(airt_time)/24 + datenum(2007,6,8,4,0,0);
air_t     = double(squeeze(ncread(ARTname,'AT_21')));
Qair_t    = double(squeeze(ncread(ARTname,'QAT_5021')));

wind_time = ncread(WINDname,'time');
wind_time = double(wind_time)/24 + datenum(2007,6,8,4,0,0);
zwind     = -double(ncread(WINDname,'depth'));
wind_spd  = double(squeeze(ncread(WINDname,'WS_401')));
Qwind_spd = double(squeeze(ncread(WINDname,'QWS_5401')));
wind_dir  = double(squeeze(ncread(WINDname,'WD_410')));
Qwind_dir = double(squeeze(ncread(WINDname,'QWD_5410')));
wind_str(Qwind_str==0) = NaN;
air_t(Qair_t==0)       = NaN;
wind_spd(Qwind_spd==0) = NaN;
wind_dir(Qwind_dir==0) = NaN;
wind_dir = wind_dir + 180; % clockwise from N [degree, from]

% interpolate meteorological measurements to the time of wave measurements
wstr = interp1(tau_time,wind_str,wave_time)';
airt = interp1(airt_time,air_t,wave_time)';
wspd = interp1(wind_time,wind_spd,wave_time)';
wdir = interp1_ang(wind_time,wind_dir,wave_time)';

% adjust wind speed to the height of half wave length
zHWL = g/4/pi./(fctr.^2);
wspdHWL = spshfttc(wspd,zwind,zHWL,airt);

F82 = Fin_from_dws(wq,fctr,wave_bw,wstr,   wdir,'P82');
F87 = Fin_from_dws(wq,fctr,wave_bw,wspdHWL,wdir,'DP87');

% split matrix into chunks due to the MATLAB size limit
nChunk = fix(ntm/10000);
jt = (0:nChunk)'*10000;
jt(1) = 1;
jt(end) = ntm + 1;

for j = 1:nChunk
    
    jl = jt(j);
    jr = jt(j+1) - 1;
                   
    [uSt(:,jl:jr),...
     vSt(:,jl:jr)] = Stokes_from_dws(wq(jl:jr,:),fctr,wave_bw,zSt);
 
    [duSt_dz(:,jl:jr),...
     dvSt_dz(:,jl:jr)] = StShear_from_dws(wq(jl:jr,:),fctr,wave_bw,ziSt);
end

%% Timetable

F82 = F82';
F87 = F87';
uSt = uSt';
vSt = vSt';
duSt_dz = duSt_dz';
dvSt_dz = dvSt_dz';

sd = timetable(datm,dahr,uSt,vSt,duSt_dz,dvSt_dz,Hs,Cp,F82,F87);

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

save([data_dir,'Papa/ocsp_wave_1hr.mat'],'SD','zSt');

%% 3-Hour average

% % bin edges and center
% tlower   = dateshift(datm(1),  'start','day');
% tupper   = dateshift(datm(end),'end',  'day');
% datm_3hr = (tlower:hours(3):tupper)';
% datm3E   = datm_3hr - hours(1) - minutes(30);
% 
% SD3 = retime(sd,datm3E,'mean');
% SD3.datm = datm_3hr;
% SD3.dahr = datenum(datm_3hr)*24; % in hours
% 
% % eliminate bins with more/less than 3 data points for 3-hour average
% SD = SD3(7:end-8,:);
% 
% save([data_dir,'Papa/ocsp_wave_3hr.mat'],'SD','zSt');

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

wave_r1 = sqrt(wave_a1n.^2 + wave_b1n.^2);
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