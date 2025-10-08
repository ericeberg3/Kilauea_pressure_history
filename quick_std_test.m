downsamplerate = 720;
load(['data/u_mm_' int2str(downsamplerate) '.mat']);
load('data/gps_time_series.mat');
load('data/GPS_seism_locations.mat', 'GPSNameList', 'GPS_llh');
load('data/hawaii_line_new.mat', 'coast_new', 'new_pit')
load('data/sdh_clean.mat')

t = decyear(sdh_clean.t);
year = floor(t);
partialYear = mod(t,1);
date0 = datenum(num2str(year),'yyyy');
date1 = datenum(num2str(year+1),'yyyy');
daysInYear = date1 - date0;
t = date0 + partialYear .* daysInYear;
t = datetime(t, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy');
% collapset = datetime(collapset, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy');

%% Figure creation
figure(1);
% Use feb 08 to Mar 15 for random walk fitting
clf;
hold on;
plot(t, sdh_clean.d(:, 1));
xline(t(1e4))
xline(t(6e4))