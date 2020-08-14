function read_ocsp_flux(cGotm_in)
%
% Read computed air-sea fluxes from the OCSP mooring
% Average the time series in one-hour and three-hours intervals
%
% cGotm_in: option to create input files for GOTM simulation
%  1 - true
%  0 - false
%
% Averaged data are saved as 'ocsp_flux_1hrPMEL.mat' and
%  'ocsp_flux_3hrPMEL.mat', respectively
%
% Note: sensible heat flux due to Rain is not included
%
% Zhihua Zheng, University of Washington, July 15 2019
% =========================================================================

%% General setting

root_dir = '~/GDrive/UW/Research/Data/';
Met_dir  = [root_dir,'Papa/Met/'];
TS_dir   = [root_dir,'Papa/TS/'];

TAUname = fullfile(Met_dir,'tau50n145w_hr.cdf');
NSWname = fullfile(Met_dir,'swnet50n145w_hr.cdf');
NLWname = fullfile(Met_dir,'lwnet50n145w_hr.cdf');
QLHname = fullfile(Met_dir,'qlat50n145w_hr.cdf');
QSHname = fullfile(Met_dir,'qsen50n145w_hr.cdf');
Ename   = fullfile(Met_dir,'evap50n145w_hr.cdf');
Pname   = fullfile(Met_dir,'rain_wspd_cor50n145w_hr.cdf');

SSTname = fullfile(TS_dir,'sst50n145w_hr.cdf');
SSSname = fullfile(TS_dir,'sss50n145w_hr.cdf');
ARTname = fullfile(Met_dir,'airt50n145w_hr.cdf');
WINname = fullfile(Met_dir,'w50n145w_hr.cdf');

%% Read variables

time  = ncread(TAUname,'time');

% wind stress
tau_x = ncread(TAUname,'TX_442');  % [N/m^2]
tau_y = ncread(TAUname,'TY_443');  % [N/m^2]
tau   = ncread(TAUname,'TAU_440'); % [N/m^2]
tauQ  = ncread(TAUname,'QTAU_5440');

% net solar radiation
nsw   = ncread(NSWname,'SWN_1495'); % [W/m^2] into the ocean
nswQ  = ncread(NSWname,'QSW_5495');

% net longwave radiation
nlw   = -ncread(NLWname,'LWN_1136'); % [W/m^2] into the ocean
nlwQ  =  ncread(NLWname,'QLW_5136');

% latent heat flux
hlb   = -ncread(QLHname,'QL_137');   % [W/m^2] into the ocean
hlbQ  =  ncread(QLHname,'QQL_5137');

% sensible heat flux
hsb   = -ncread(QSHname,'QS_138');   % [W/m^2] into the ocean
hsbQ  =  ncread(QSHname,'QQS_5138');

% evaporation
evap  = ncread(Ename,'E_250');  % [mm/hr]
evapQ = ncread(Ename,'QE_5250');

% rain
rain  = ncread(Pname,'RN_485'); % [mm/hr]
rainQ = ncread(Pname,'QRN_5485');
train = ncread(Pname,'time');

% sea surface temperature
sst  = ncread(SSTname,'T_25'); % [C]
sstQ = ncread(SSTname,'QT_5025');
tsst = ncread(SSTname,'time');

% sea surface salinity
sss  = ncread(SSSname,'S_41'); % [PSU]
sssQ = ncread(SSSname,'QS_5041');
tsss = ncread(SSSname,'time');

% air temperature
art  = ncread(ARTname,'AT_21'); % [C]
artQ = ncread(ARTname,'QAT_5021');
tart = ncread(ARTname,'time');

% surface wind
uwin =  ncread(WINname,'WU_422');
vwin =  ncread(WINname,'WV_423');
zwin = -ncread(WINname,'depth'); % [m]
winQ =  ncread(WINname,'QWS_5401');
twin =  ncread(WINname,'time');

%% double precisions

time  = double(time);
train = double(train);
tsst  = double(tsst);
tsss  = double(tsss);
tart  = double(tart);
twin  = double(twin);

tau_x = double(squeeze(tau_x));
tau_y = double(squeeze(tau_y));
tau   = double(squeeze(tau));
nsw   = double(squeeze(nsw));
nlw   = double(squeeze(nlw));
hlb   = double(squeeze(hlb));
hsb   = double(squeeze(hsb));
evap  = double(squeeze(evap));
rain  = double(squeeze(rain));
sst   = double(squeeze(sst));
sss   = double(squeeze(sss));
art   = double(squeeze(art));
uwin  = double(squeeze(uwin));
vwin  = double(squeeze(vwin));

tauQ  = double(squeeze(tauQ));
nswQ  = double(squeeze(nswQ));
nlwQ  = double(squeeze(nlwQ));
hlbQ  = double(squeeze(hlbQ));
hsbQ  = double(squeeze(hsbQ));
evapQ = double(squeeze(evapQ));
rainQ = double(squeeze(rainQ));
sstQ  = double(squeeze(sstQ));
sssQ  = double(squeeze(sssQ));
artQ  = double(squeeze(artQ));
winQ  = double(squeeze(winQ));

%% Select high quality data

% quality codes have category 0 and 1
tau_x(tauQ==0) = NaN;
tau_y(tauQ==0) = NaN;
tau(tauQ==0)   = NaN;
nsw(nswQ==0)   = NaN;
nlw(nlwQ==0)   = NaN;
hlb(hlbQ==0)   = NaN;
hsb(hsbQ==0)   = NaN;
evap(evapQ==0) = NaN;
rain(rainQ==0) = NaN;
uwin(winQ==0)  = NaN;
vwin(winQ==0)  = NaN;
art(artQ==0)   = NaN;

% quality codes have category 0, 1, 2, 3, 5
sst(sstQ==0 | sstQ==5) = NaN;
sss(sssQ==0 | sssQ==5) = NaN;

%% Adjust timestamp

time  = datenum(2007,6,8,5,0,0) + time/24;
train = datenum(2007,6,8,5,0,0) + train/24;
tsst  = datenum(2007,6,8,0,0,0) + tsst/24;
tsss  = datenum(2007,6,8,0,0,0) + tsss/24;
tart  = datenum(2007,6,8,0,0,0) + tart/24;
twin  = datenum(2007,6,8,4,0,0) + twin/24;
twin  = round(twin*24*3600)/24/3600; % round to nerest second

dahr  = 24*time; % in hours
datm  = datetime(time, 'ConvertFrom','datenum');
drain = datetime(train,'ConvertFrom','datenum');
dsst  = datetime(tsst, 'ConvertFrom','datenum');
dsss  = datetime(tsss, 'ConvertFrom','datenum');
dart  = datetime(tart, 'ConvertFrom','datenum');
dwin  = datetime(twin, 'ConvertFrom','datenum');

% truncate other time series to the time of wind stress
sInx = find(dsst == datm(1));
eInx = find(dsst == datm(end));
sst  = sst(sInx:eInx);

sInx = find(dsss == datm(1));
eInx = find(dsss == datm(end));
sss  = sss(sInx:eInx);

sInx = find(drain == datm(1));
rain = rain(sInx:end);
rain = [rain; nan(length(tau)-length(rain),1)];

sInx = find(dart == datm(1));
eInx = find(dart == datm(end));
art  = art(sInx:eInx);

sInx = find(dwin == datm(1));
eInx = find(dwin == datm(end));
uwin = uwin(sInx:eInx);
vwin = vwin(sInx:eInx);
wspd = sqrt(uwin.^2 + vwin.^2);
u10  = spshfttc(wspd,zwin,10,art);

SF = timetable(datm,dahr,...
     tau_x,tau_y,tau,u10,nsw,nlw,hlb,hsb,evap,rain,sst,sss);

save([root_dir,'Papa/ocsp_flux_1hrPMEL.mat'],'SF');

%% GOTM input files

if cGotm_in
    
gotmdata_root = '~/Documents/GitHub/GOTM/gotmwork/data/';
basecase      = [gotmdata_root,'OCSPapa_20070608-20190616/'];

TAUname  = [basecase,'tau_file.dat'];
HFname   = [basecase,'heatflux_file.dat'];
HFSWname = [basecase,'heatflux_swr.dat']; % heat flux including radiation
RAINname = [basecase,'precip_file.dat'];
SSSname  = [basecase,'sss_file.dat'];
SSTname  = [basecase,'sst_file.dat'];
SWRname  = [basecase,'swr_file.dat'];

hf       = hlb + hsb + nlw;
hfsw     = hf  + nsw;
tau_xy   = [tau_x tau_y];
rainG    = rain/1000/3600; % GOTM needs precipitation in [m/s] !
datm_str = string(datestr(datm,'yyyy-mm-dd HH:MM:SS'));

write_gotm_flux(TAUname, tau_xy(~isnan(tau),:),datm_str(~isnan(tau)))
write_gotm_flux(HFname,  hf(~isnan(hf)),       datm_str(~isnan(hf)))
write_gotm_flux(HFSWname,hfsw(~isnan(hfsw)),   datm_str(~isnan(hfsw)))
write_gotm_flux(RAINname,rainG(~isnan(rain)),  datm_str(~isnan(rain)))
write_gotm_flux(SSSname, sss(~isnan(sss)),     datm_str(~isnan(sss)))
write_gotm_flux(SSTname, sst(~isnan(sst)),     datm_str(~isnan(sst)))
write_gotm_flux(SWRname, nsw(~isnan(nsw)),     datm_str(~isnan(nsw)))
end

end