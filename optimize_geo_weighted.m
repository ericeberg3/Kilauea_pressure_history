%% This script uses green's functions created from an optimized HMM and SC to invert for the pressure history

% Loading in data and setting the downsample rate (from 5s to hourly)
downsamplerate = 720;
load(['data/u_mm_full_duration_' int2str(downsamplerate) '_interp.mat']);
load('data/gps_time_series.mat');
load('data/GPS_seism_locations.mat', 'GPSNameList', 'GPS_llh');
load('data/hawaii_line_new.mat', 'coast_new', 'new_pit')
load('data/sdh_clean.mat')

%% Creating GPS names and setting up x and y coordinates
GPSNameList = ["69FL","92YN","AHUP","BDPK","BYRL","CNPK","CRIM","DEVL","OUTL","PUHI","PWRL","UWEV","V120","VSAS", "CALS"];
tiltNameList = ['SDH'];

t = downsample(universalt, downsamplerate);
ntime = length(t);

% convert u from mm to m
u = u .* 1e-3;

% Set outliers to NaN
for i = 1:size(u, 1)
    u(i, :, :) = filloutliers(squeeze(u(i, :, :)), NaN, "movmedian", 300, 2);
end 

% Center the displacements at the beginning of each time series at 0. Done 
% by taking the median of the first 20 non-nan data points for each station
u = CenterDisps(20, u, GPSNameList, ntime);

clearvars interpu unfilledu nanu meansize rtnet_ka windowsize startcind new_pit endcind disps;

nstat = length(GPSNameList); % Number of GPS stations

% station x and y positions
xy = llh2local(GPS_llh, [-155.2784, 19.4073]);
x = xy(1,:) * 1000;
y = xy(2,:) * 1000;
z = zeros(1, nstat);

%% Align tilt data to displacment data
tilts = interp1(decyear(sdh_clean.t), sdh_clean.d, universalt);
tilts = medfilt1(tilts, downsamplerate/40, [], 1, 'omitnan', 'truncate');
tilts = interp1(universalt, tilts, t, 'linear');

% Tilt standard dev. is calculated by taking the mean of the E and N
% components from April 11th - April 18th 2018
tilts = (tilts - tilts(1, :));
tilts(1, :) = 0;
tiltstd = std(sdh_clean.d(1e5:1.1e5, 1:2), 1);
tiltstd = mean(tiltstd);

%% Manually set theta offset
% dtheta = -deg2rad(0);
% % dphi = tiltoffsets(2);
% 
% disp("Theta offset: " + dtheta)
% disp("phi offset: " + dphi)


%% Cleaning up stations with too many nan values

% Removing CALS
u(end, :, :) = [];
GPSNameList(end) = [];
nstat = nstat - 1;

%% Setting station locations and eliminating faulty stations 

% Only get the displacement data from the first and last timestamp for
% use in chamber optimization. u1d & tiltreduced are the long term
% displacements / tilts
finalindex = 114;
u1d = u(:, :, end-finalindex);

nanstatend = isnan(u1d(:, 2));
nanstatbeginning = isnan(u(: , 1, 1));

%%% Uncomment if we don't want to include CRIM, UWEV, BYRL %%% 
nanstatend(find(contains(GPSNameList,'CRIM')), :, :) = 1;
nanstatend(find(contains(GPSNameList,'UWEV')), :, :) = 1;
nanstatend(find(contains(GPSNameList,'BYRL')), :, :) = 1;

% Delete nan stations from the data
u1d = u1d(~nanstatend, :);
tiltreduced = tilts(end-finalindex, :) - tilts(1, :);
xopt = x(~nanstatend);
yopt = y(~nanstatend);
zopt = zeros(1, length(xopt));
nanstatbeginning = nanstatbeginning(~nanstatend);
offsets = 0.5 .* ones(1, sum(nanstatbeginning == 1));

% Splitting overall u into more readable ux, uy, uz
ux = squeeze(u(:, 1, :));
uy = squeeze(u(:, 2, :));
uz = squeeze(u(:, 3, :));

%% Setting the tilt location and NPIT location for optimization

% Tilt location
xyll = [-155.2937; 19.3905];
xy = llh2local(xyll, [-155.2784, 19.4073]) * 1000;
xtilt = xy(1);
ytilt = xy(2);

% NPIT location
npitloc = coord('NPIT', 'llh');
npitloc = llh2local(npitloc(1:2), [-155.2784, 19.4073]) * 1000;
vertscale = 4e3;
radscale = 4e3;

% Setting up weights for optimization
invStdPWRL = 1./std(squeeze(u(11, :, :)), 0, 2, "omitmissing");
invStdPWRL = invStdPWRL(:);

%% Optimizing SC geometry
mSCguess = [277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, 1e7];


% Use taiyi's priors except for strike angle and delta p
taiyi_parameters = [1600.79, 914.47, 90, 0, 70.3, 183, -1940, -3e6, ... 
    277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, -10e6];
% taiyi_parameters(1:2) = [2015, 1151]; % 2x volume
taiyi_parameters(7) = -2400; % 2x volume

% taiyi_parameters(1:2) = [2540, 1450]; % 4x volume
% taiyi_parameters(7) = -2880; % 4x volume
% taiyi_parameters_flat_SC = [1600.79, 914.47, 90, 0, 70.3, 183, -1940, -3e6, ... 
    % 277.01, 1621.47, 90, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, -10e6];
% lb = [-4e6, 60, 126, mSCguess(5) - 150, mSCguess(6) - 200, -3.85e3, -2e7];
% ub = [-2e6, 120, 146, mSCguess(5) + 150, mSCguess(6) + 200, -3.43e3, -3e6];

% Let the parameters be mostly free
% dpHMM, vertical semi-diameter, horizontal semi-diameter, dip, strike, x1,
% x2, x3, dpSC
% lb = [-5e6, 100, 1e3, 60, 126, mSCguess(5) - 2e3, mSCguess(6) - 2e3, -4.5e3, -2e7];
% ub = [0, 1e3, 1e4, 90, 146, mSCguess(5) + 2e3, mSCguess(6) + 2e3, -2.5e3, 0];
lb = [-10e6, 100, 1e3, 62, -1e7];
ub = [0, 1e3, 1e4, 90, -9e6];

% optParams = optimize_SC_MCMC(mSCguess, taiyi_parameters, lb, ub, xopt, xtilt, yopt, ytilt, zopt, u1d, invStdPWRL, tiltstd, tiltreduced, nanstatbeginning);
% Old optimization routine
optParams = optimize_SC_bayes(mSCguess, taiyi_parameters, lb, ub, xopt, xtilt, yopt, ytilt, zopt, u1d, invStdPWRL, tiltstd, tiltreduced, nanstatbeginning);
optimizedM = [taiyi_parameters(1:7), optParams.dpHMM, optParams.vert_sd, optParams.horiz_sd, optParams.dip, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, optParams.dpSC];
% optimizedM =[1600.79000000000	914.470000000000	90	0	70.3000000000000	183	-1940	-2999999.99992412	277.010000000000	1621.47000000000	89.9999774942710	136	162.745313651300	-479.303502915132	-4620.30457577829   -12317057.4035003];

% optimizedM = taiyi_parameters_flat_SC;
offsets = optimizedM(17:end);
optimizedM = optimizedM(1:16);


%% Now optimize for the tilt orientation:
% f = @(t)green_residuals_tilt(t, [optimizedM, offsets], [xopt, xtilt], [yopt, ytilt], [zopt, 0], u1d, ...
%     [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3),1/((tiltstd)), 1/((tiltstd))], tiltreduced(1:2), nanstatbeginning);
% problem = createOptimProblem('fmincon','objective', f,'x0', 0,'lb',-pi/4, ...
%     'ub',pi/4,'options',options);
% ms = MultiStart;
% [dtheta, res] = run(ms,problem,20);
dtheta = deg2rad(0);
disp("Optimized theta = " + rad2deg(dtheta))

% tilts(:, 1) = 1e6 * atan(tan(1e-6 * tilts(:, 1)) * cos(dtheta) - tan(1e-6 * tilts(:, 2)) * sin(dtheta));
% tilts(:, 2) = 1e6 * atan(tan(1e-6 * tilts(:, 2)) * sin(dtheta) + tan(1e-6 * tilts(:, 2)) * cos(dtheta));

% Adding the optimized offsets to the GPS data
j = 1;
for i = 1:length(nanstatbeginning)
    if(nanstatbeginning(i) == 1)
        ux(i, :) = ux(i, :); % + offsets(j);
        uy(i, :) = uy(i, :); % + offsets(j);
        uz(i, :) = uz(i, :); % + offsets(j);
        j = j + 1;
    end
end

clearvars i j xopt yopt zopt options lb ub resnorm resdpos dsize dangle A b res;

%% Creating greens functions using the new parameters for M
% dtheta = -0.9 * pi/4; % In radians the offset angle of the tiltmeter
% dphi = 3 * pi/4;
[gHMM, gSC] = creategreens(optimizedM(1:8), optimizedM(9:end));
[gTiltHMM, gTiltSC] = createtiltgreens(optimizedM(1:8), optimizedM(9:end), dtheta, false);

%% Calculate the pressure at each time step
% Setting up green's functions and formatting displacements/tilts properly
gHMMflat = gHMM';
gHMMflat = gHMMflat(:);
gSCflat = gSC';
gSCflat = gSCflat(:);
gTiltHMMflat = gTiltHMM';

ntime = max(size(u(1,1,:)));

tiltx = tilts(:, 1); % * (cos(dtheta) - sin(dtheta))/cos(dphi); 
tilty = tilts(:, 2); % * (cos(dtheta) + sin(dtheta))/cos(dphi);

i = 1;
dispstd = [invStdPWRL(1) .* ones(1, length(ux(:, i))), invStdPWRL(2) .* ones(1, length(ux(:, i))), invStdPWRL(3) .* ones(1, length(ux(:, i))), 1/((tiltstd)), 1/((tiltstd))];
temp = TimeDependentLSQtilt(gHMMflat, gSCflat, gTiltHMM, gTiltSC, ux, uy, uz, tiltx, tilty, dispstd, GPSNameList);
dp = real([temp(1:length(tiltx)); temp(length(tiltx) + 1:end)])';
clear temp i;

% Fill outliers, which occur due to gaps in data throwing off LSQ suddenly
[dp(:, 2), TF] = filloutliers(dp(:, 2), "makima", "movmedian", 100, 1);
dp(TF, 1) = NaN;
dp(:, 1) = fillmissing(dp(:, 1), "makima", 1);

%% Creating a matrix to store all predicted displacements + tilts in var called usim
usim = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));

for i = 1:max(size(dp))
    % **** SETTING THE SIMULATED DISPLACEMENTS ****
    usim(i, :, 1:nstat) =  (gHMM .* dp(i, 1)) + (gSC .* dp(i, 2)) ;
    usim(i, :, nstat + 1) = [gTiltHMM .* dp(i, 1) + gTiltSC .* dp(i, 2), 0];
end

%% Extract collapse amplitudes
collapset = D.gps_info.t_events(26:end);
collapset = decyear(datetime(collapset, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy'));
% [ampHMM, ampSC] = ExtractCollapseAmplitude([dp(:, 1)'; dp(:, 2)'], t, collapset, t(3) - t(2));

%%
makeplots(x, y, z, u, ux, uy, uz, tiltx, tilty, usim, t, finalindex, collapset, dp, optimizedM, ...
    GPSNameList, gTiltHMM, gTiltSC, xtilt, ytilt, tiltreduced, radscale, mSCguess, coast_new, dtheta, 3);
