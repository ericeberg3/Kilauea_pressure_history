%% This script uses green's functions created from an optimized HMM and SC to invert for the pressure history

% Loading in data and setting the downsample rate (from 5s to hourly)
downsamplerate = 720;
% u_mm_full_duration was created from the script generatecontinuousdisplacements in the folder time_series 
% load(['data/u_mm_full_duration_' int2str(downsamplerate) '_interp.mat']);
load(['data/u_mm_full_duration_' int2str(downsamplerate) '_interp.mat']);
load('data/gps_time_series.mat');
load('data/GPS_seism_locations.mat', 'GPSNameList', 'GPS_llh');
load('data/hawaii_line_new.mat', 'coast_new', 'new_pit')
load('data/sdh_clean.mat')
daily_GPS = readtable('data/daily_gps_mult_stations.csv');

addpath('./Plotting functions/')

%% Creating GPS names and setting up x and y coordinates
GPSNameList = ["69FL","92YN","AHUP","BDPK","BYRL","CNPK","CRIM","DEVL","OUTL","PUHI","PWRL","UWEV","V120","VSAS", "CALS"];
tiltNameList = ['SDH'];

% Downsampling the time vector
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
u = CenterDisps(1, 20, u, GPSNameList, ntime);

clearvars interpu unfilledu nanu meansize rtnet_ka windowsize startcind new_pit endcind disps;

nstat = length(GPSNameList); % Number of GPS stations

% station x and y positions
xy = llh2local(GPS_llh, [-155.2784, 19.4073]);
x = xy(1,:) * 1000;
y = xy(2,:) * 1000;
z = zeros(1, nstat - 1);

% Removing CALS (not used in this study)
u(end, :, :) = [];
GPSNameList(end) = [];
nstat = nstat - 1;

%% Set up daily GPS data for geometry inversion saved in the u1d array
nanstatend = false(1, length(GPSNameList));

% This time interval gives us the most GPS data
t0 = datetime(2018,5,19,'TimeZone','UTC');
t1 = datetime(2018,8,10,'TimeZone','UTC');

% Get the start and end index of from the universal time vector
gpsDates = datetime(daily_GPS.Date, ...
                    'InputFormat',"MM/dd/yyyy HH:mm", ...
                    'TimeZone',"UTC") + calyears(2000);
[~, startind] = min(abs(gpsDates - t0));
[~, finalind] = min(abs(gpsDates - t1));

% Get daily GPS Data from PWRL as a proxy for white noise
PWRLdata = [daily_GPS{startind:finalind, "PWRL_east"}, daily_GPS{startind:finalind, "PWRL_north"}, daily_GPS{startind:finalind, "PWRL_up"}];
daily_inv_std = 1./std(PWRLdata, "omitnan");

% Construct the long term displacement array (u1d = 3xNstat array)
u1d = zeros(length(GPSNameList), 3);
for i = 1:length(GPSNameList)
    station = char(GPSNameList(i));
    % Need to check if the station starts with a number, if so the
    % name list places an 'x' in front of the name
    check = all(ismember(station(1), '0123456789'));

    % Excluding CRIM, UWEV, BYRL from subsequent inversions by flagging
    % them in the nanstatend array
    if(GPSNameList(i) == "CRIM" || GPSNameList(i) == "UWEV" || GPSNameList(i) == "BYRL"); nanstatend(i) = true; end
    try
        if(check == 0)
            % We can populate the u1d array using the name found
            u1d(i, 1) = daily_GPS{finalind, GPSNameList(i) + "_east"} - daily_GPS{startind, GPSNameList(i) + "_east"};
            u1d(i, 2) = daily_GPS{finalind, GPSNameList(i) + "_north"} - daily_GPS{startind, GPSNameList(i) + "_north"};
            u1d(i, 3) = daily_GPS{finalind, GPSNameList(i) + "_up"} - daily_GPS{startind, GPSNameList(i) + "_up"};
        else
            % If it's a number we must add an 'x' in front
            u1d(i, 1) = daily_GPS{finalind, "x" + GPSNameList(i) + "_east"} - daily_GPS{startind, "x" + GPSNameList(i) + "_east"};
            u1d(i, 2) = daily_GPS{finalind, "x" + GPSNameList(i) + "_north"} - daily_GPS{startind, "x" + GPSNameList(i) + "_north"};
            u1d(i, 3) = daily_GPS{finalind, "x" + GPSNameList(i) + "_up"} - daily_GPS{startind, "x" + GPSNameList(i) + "_up"};
        end
    catch
        % If we have errors getting data just mark the station as nan
        nanstatend(i) = true;
    end
    % Print out the stations that are nan at the beginning
    if(isnan(u1d(i, 1))); nanstatend(i) = true; disp(station); end
end
% Delete all irrelevant stations. Can change this line if we want to
% include more stations
u1d = u1d(~nanstatend, :);

clear PWRLdata daily_GPS

%% Get random walk noise value from data
rw_stddev = GetRandomWalk(sdh_clean);

%% Align tilt data to displacment data
tilts = interp1(decyear(sdh_clean.t), sdh_clean.d, universalt);
tilts = interp1(universalt, tilts, t, 'linear');

% Tilt standard dev. is calculated by taking the mean of the E and N
% components from April 11th - April 18th 2018
tilts = (tilts - tilts(1, :));
tilts(1, :) = 0;
% tiltstd is just used for chi2 stat. Not really a good metric.
tiltstd = std(sdh_clean.d(1e5:1.1e5, 1:2), 1);
tiltstd = mean(tiltstd);

clear universalt

%% Setting station locations and eliminating faulty stations 

% Only get the displacement data from the first and last timestamp for
% use in chamber optimization. u1d & tiltreduced are the long term
% displacements / tilts
finalindex = 149; % End of the co-collapse sequence

% Start tilts at 0 and cut out unused timesteps
tiltreduced = tilts(end-finalindex, :) - tilts(1, :);

% Delete nan stations from the data
xopt = x(~nanstatend);
yopt = y(~nanstatend);
zopt = zeros(1, length(xopt));

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
vertscale = 5e3;
radscale = 5e3;

% Setting up GPS stddev for LSQ later
invStdPWRL = 1./std(squeeze(u(11, :, :)), 0, 2, "omitmissing");
invStdPWRL = invStdPWRL(:);

%% Optimizing geometry
% Setting up parameters from Wang et al. 2021 for reference
% Use taiyi + kyle for HMM 
% taiyi_parameters = [1600.79, 914.47, 90, 0, 50, 200, -2.18e3, -4e7, ... 
%      277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, -10e6];
taiyi_parameters = [1600.79, 914.47, 90, 0, 50, 200, -2.18e3, -4.5e6, ... 
     277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, -1.5e7];

% Setting up parameters from Roman et al. 2021 for reference
aspect_ratio_HMM = 1;
aspect_ratio_SC = 1;
roman_HMM_vol = 9.1 * 1e9; 
roman_SC_vol = 10.5 * 1e9;
vert_sd_HMM = (3/(4*pi) * roman_HMM_vol * (aspect_ratio_HMM^2))^(1/3);
vert_sd_SC = (3/(4*pi) * roman_SC_vol * (aspect_ratio_SC^2))^(1/3);
roman_parameters = [vert_sd_HMM, vert_sd_HMM/(aspect_ratio_HMM), 90, 0, 0, 0, -900 - vert_sd_HMM, -3e6, ...
    vert_sd_SC, vert_sd_SC/(aspect_ratio_SC), 90, 136, -200, -1.9e3,-3.2e3 - vert_sd_SC, -1e6];

% Variable sequence: FIX VOLUME BOUNDS
% ["dpHMM_insar", "dpHMM_gps", "volHMM", 'xHMM', 'yHMM', "zHMM" (from top of spheroid), alpha_HMM
% xSC, ySC, dSC, "SC aspect ratio", "dip", strike, "dpSC_insar", "dpSC_gps", volSC]
% yHMM lb used to be -0.5e3 + npitloc(2)
% Old:
% lb = [-5e7, -5e7, 1e8, -5e2, -0.55e3 + npitloc(2), -2e3, 0.8, ...
%     0.7e3, -2.9e3, -4.6e3, 0.05, 40, 100, -5e7, -5e7, 2.0e9]; 
% ub = [1e6, 1e6, 3e10, 0.5e3 + npitloc(1), 0.5e3 + npitloc(2), -3e2, 2.0, ...
%      1.8e3, -2.0e3, -2.9e3, 1, 90, 180, 1e6, 1e6, 20e9]; 

% Inferred from kyle / taiyi (see spreadsheet):
% lb = [-5e7, -5e7, 1e8, -100, 0, -16e2, 0.8, ...
%     -2.7e3, -2.8e3, -4.7e3, 0.1, 45, 0, -5e7, -5e7, 2.0e9]; 
% ub = [1e6, 1e6, 3e10, 227, 1050, -7.5e2, 1.8, ...
%      2.2e3, -0.45e3, -2.7e3, 1, 90, 180, 1e6, 1e6, 20e9]; 


lb = [-1e8, -1e8, -100, 0, -16e2, 0.8, ...
    -2.7e3, -2.8e3, -4.7e3, 0.1, 45, 0, -1e8, -1e8]; 
ub = [1e8, 1e8, 227, 1050, -7.5e2, 1.8, ...
     2.2e3, -0.45e3, -2.7e3, 1, 90, 180, 1e8, 1e8]; 

% Save figures for export
saveFigs = true;

%% Loading in InSAR data for geometry inversion
% All insar data was created using scripts in the "InSAR processing" folder
% see the readme for more detail on how the data was created

% Displacement data
insar_data_ascending = readmatrix('data/ascending.txt', 'Delimiter', '\t');
insar_data_descending = readmatrix('data/descending.txt', 'Delimiter', '\t');

% Covariance data
cov_asc = load("Data/asc_cov.mat", "insarcov_quad").insarcov_quad;
cov_desc = load("Data/desc_cov.mat", "insarcov_quad").insarcov_quad;

% Look vectors
load 'data/look_asc.mat';
load 'data/look_desc.mat';

% local x, y bounds for insar data
insarbnd = [-1e4, -1e4; 1e4, 1e4];
look = [avglook_asc, avglook_desc];

% Process asc track
insarxy=llh2local(insar_data_ascending(:,1:2)', [-155.2784, 19.4073])'.*1e3;
in_bounds = (insarxy(:,1) >= insarbnd(1,1)) & (insarxy(:,1) <= insarbnd(2,1)) & ...
            (insarxy(:,2) >= insarbnd(1,2)) & (insarxy(:,2) <= insarbnd(2,2));
% X, Y, coordinates of each insar data point
insarx_asc = insarxy(in_bounds,1)';
insary_asc = insarxy(in_bounds,2)';
% Cropped displacement data points
insaru_asc = -insar_data_ascending(in_bounds, 3)'; %negative bc this direction is towards the satellite
% Blocks represents the size of each block (we have already downsampled)
blocks_asc = insar_data_ascending(in_bounds, 4)';
% Preemptively taking the inverse of the covariance matrix
cov_asc = cov_asc(in_bounds,in_bounds);
cinv_asc = pinv(cov_asc);

% Process descending track
insarxy=llh2local(insar_data_descending(:,1:2)', [-155.2784, 19.4073])'.*1e3;
in_bounds = (insarxy(:,1) >= insarbnd(1,1)) & (insarxy(:,1) <= insarbnd(2,1)) & ...
            (insarxy(:,2) >= insarbnd(1,2)) & (insarxy(:,2) <= insarbnd(2,2));
% X, Y, coordinates of each insar data point
insarx_desc = insarxy(in_bounds,1)';
insary_desc = insarxy(in_bounds,2)';
% Cropped displacement data points
insaru_desc = -insar_data_descending(in_bounds, 3)'; %negative bc this direction is towards the satellite
% Blocks represents the size of each block (we have already downsampled)
blocks_desc = insar_data_descending(in_bounds, 4)';
% Preemptively taking the inverse of the covariance matrix
cov_desc = cov_desc(in_bounds,in_bounds);
cinv_desc = pinv(cov_desc);

% Combining asc and desc track inverse cov matrix for the geometry inversion
cinv_full = blkdiag(cinv_asc, cinv_desc);

% Collect full insar data into new variables
insar_lengths = [length(insaru_asc), length(insaru_desc)];
insaru_full = [insaru_asc, insaru_desc];
insarx = [insarx_asc, insarx_desc];
insary = [insary_asc, insary_desc];
block_size = [blocks_asc, blocks_desc];

%% MCMC Static inversion #1 - Volume change
% paramNames = {'dpHMM_insar', 'dpHMM_gps', 'volHMM', 'xHMM', 'yHMM', 'dHMM', 'alphaHMM' ...
%     'xSC', 'ySC', 'dSC', 'alphaSC', 'dipSC', 'strikeSC', 'dpSC_insar', 'dpSC_gps', 'volSC'};
paramNames = {'dvHMM_insar', 'dvHMM_gps', 'xHMM', 'yHMM', 'dHMM', 'alphaHMM' ...
    'xSC', 'ySC', 'dSC', 'alphaSC', 'dipSC', 'strikeSC', 'dvSC_insar', 'dvSC_gps'};
ntrials = 1e5; % Customize to get convergence

% Testing GPS and prior weights.
% gps_weights = linspace(4e1, 8e1, 10);
% prior_weights = linspace(3e2, 1e3, 10);
gps_weights = linspace(8e1, 6e2, 15);
prior_weights = linspace(7e1, 1e5, 10);

gps_weight = 6.7e1; % Optimal weight based on L curve
prior_weight = 0;%5.3e2; % Optimal weight based on L curve
% burn = 0.5e3;
burn = 4e3;

delete(gcp('nocreate'));

% --- Set if we want to run MCMC, plot the L curve, and the type of L curve
% to create --- %
runMCMC = false;
run_L_curve = false;
l_curve_type = "gps"; % 'prior' to test prior weights, 'gps' to test gps weights
% for l_curve_type = ["prior", "gps"]
% For loop to create L curve
n_l_curve = length(gps_weights);
if(~run_L_curve)
    n_l_curve = 1; 
else
    % Set up log likelihood lists for L curve creation. 
    gps_l2s = zeros(1,length(gps_weights));
    insar_l2s = zeros(1,length(gps_weights));
    posteriors_list = zeros(length(paramNames), ntrials - burn + 1, length(gps_weights));
    optParams_list = zeros(length(paramNames), length(gps_weights));
    l_curve_points = zeros(3,length(gps_weights));
end
start_params = taiyi_parameters;
% Move taiyi depth down a bit
start_params(7) = start_params(7) - 300;
for i = 1:n_l_curve
    prior_params = start_params;
    
    % Check if we are doing L curve analysis. If so set weight
    % appropriately
    if(run_L_curve)
        if(l_curve_type == "prior"); prior_weight = prior_weights(i); gps_weight = 6.7e1;
        else; prior_weight = 2.1e2; gps_weight = gps_weights(i); end
    end
        
    if(runMCMC)
        [optParams, posterior, L_keep, gps_l2, insar_l2, prior_l2] = optimize_SC_MCMC(prior_params, lb, ub, xopt, ...
            yopt, zopt, u1d', insarx, insary, insaru_full', look, insar_lengths, sparse(cinv_full), daily_inv_std, ...
            nanstatend, ntrials, gps_weight, prior_weight, paramNames, burn, false, saveFigs); % subsample set to true
        start_params = get_full_m(taiyi_parameters, real(optParams'), true, "insar");

        if(~run_L_curve)
            % save Data/MCMC_1e6_SCvol_topbnd.mat optParams posterior L_keep gps_l2 insar_l2 prior_l2;
            save Data/MCMC_1e6_dVinversion_noprior.mat optParams posterior L_keep gps_l2 insar_l2 prior_l2;
        else
            l_curve_points(1,i) = gps_l2;
            l_curve_points(2,i) = insar_l2;
            l_curve_points(3,i) = prior_l2;
    
            posteriors_list(:, :, i) = posterior;
            optParams_list(:, i) = optParams;
        end
    else
        % load Data/MCMC_1e6_SCvol_topbnd.mat;
        load Data/MCMC_1e6_dVinversion.mat;
        [~, ind] = max(L_keep);
        optParams = posterior(:, ind)';
    end
    
end
% Save and plot l curve data
if(run_L_curve)
    save("Data/l_curve_data_" + l_curve_type + "_prior_210.mat", "l_curve_points", "l_curve_type", "prior_weights", "gps_weights");
    % load Data/l_curve_data_prior.mat;
    posterior = squeeze(posteriors_list(:, :, end));
    optParams = optParams_list(:, end)';
    disp(gps_weights(end));
    plotLcurve(l_curve_points, l_curve_type, prior_weights, gps_weights);
end

% end

% Get the full geometry parameters based on the optimization results:
disp("GPS L2: " + gps_l2 + " InSAR L2: " + insar_l2);
optimizedM = get_full_m(taiyi_parameters, optParams, true, "insar");


%%
% Show parameter table
disp(array2table(optParams(:).', 'VariableNames',paramNames));
% Make sure parameters are real
optParams = real(optParams);

%% MCMC inversion #2: Volume + pressure change inversion
paramNames2 = {'volHMM', 'volSC', 'dvHMM_insar', 'dvHMM_gps', 'dvSC_insar', 'dvSC_gps'};

% Update priors from run 1 posterior:
load Data/paramDists.mat paramDists;
vars_to_transfer = {'dvHMM_insar', 'dvHMM_gps', 'dvSC_insar', 'dvSC_gps'};

figure(1); clf;
for i = 1:length(vars_to_transfer)
    pName = vars_to_transfer{i};

    % Find the index of this parameter in the Run 1 paramNames list
    % (We assume paramNames from Run 1 exists in workspace)
    idx = find(strcmp(paramNames, pName));
    if isempty(idx)
        warning('Parameter %s not found in Run 1 results.', pName);
        continue;
    end

    % Extract the samples
    samples = posterior(idx, :);
    
    % 3. Fit a Kernel Distribution (Smooths the histogram into a pdf)
    % This creates a probability object that supports the .pdf() method
    pd = fitdist(samples', 'Normal'); 
    subplot(2,2,i);
    histogram(samples, 50, 'Normalization','pdf');
    hold on;
    x_grid = linspace(min(samples), max(samples), 100); % 100 points for smoothness
    y_vals = pdf(pd, x_grid);
    plot(x_grid, y_vals, 'LineWidth', 2, 'Color', 'r');

    % 4. Update the paramDists structure
    paramDists.(pName).dist = pd;
    paramDists.(pName).samples = samples; 
    paramDists.(pName).family = 'prob.KernelDistribution';
    
    % Optional: Print check
    fprintf('  -> %s updated (Mean: %.2e)\n', pName, mean(samples));
end


save Data/paramDists.mat paramDists;
%%
% Run second MCMC
lb2 = [1e8, 2e9, -1e8, -1e8, -1e8, -1e8];
ub2 = [3e10, 2e10, 0, 0, 0, 0];
ntrials2 = 1e6;
runMCMC = false;
% x_guess = [4e9, 3e9, optParams(1), optParams(2), optParams(13), optParams(14)];
if(runMCMC)
    [optParams2, posterior2, L_keep2, ~, ~, ~] = optimize_SC_MCMC(optimizedM, lb2, ub2, xopt, ...
                yopt, zopt, u1d', insarx, insary, insaru_full', look, insar_lengths, sparse(cinv_full), daily_inv_std, ...
                nanstatend, ntrials2, gps_weight, 5.3e2, paramNames2, burn, true, saveFigs); % subsample set to true
    save Data/MCMC2_1e6_vol_inversion.mat optParams2 posterior2 L_keep2;
else
    load Data/MCMC2_1e6_vol_inversion.mat
end
% Adjust optimizedM to reflect the new calc volume + volume change
HMM_AR = optimizedM(1)/optimizedM(2);
HMM_volume = optParams2(1);
vert_sd = (3/(4*pi) * HMM_volume * (HMM_AR^2))^(1/3);
horiz_sd = vert_sd/(HMM_AR);
optimizedM(1:2) = [vert_sd, horiz_sd];

SC_AR = optimizedM(9)/optimizedM(10);
SC_volume = optParams2(2);
vert_sd = (3/(4*pi) * SC_volume * (SC_AR^2))^(1/3);
horiz_sd = vert_sd/(SC_AR);
optimizedM(9:10) = [vert_sd, horiz_sd];

optimizedM(8) = optParams2(3); optimizedM(16) = optParams2(5);

%% New histogram plotting
paramNames3 = {'dpHMM_insar', 'dpHMM_gps', 'dpSC_insar', 'dpSC_gps'};
lb3 = [-5e7, -5e7, -5e7, -5e7];
ub3 = [0, 0, 0, 0];
% 1. Pre-allocate pressure arrays (N_samples x 4: HMM_insar, HMM_gps, SC_insar, SC_gps)
% We assume pFromV returns [dpHMM_insar, dpHMM_gps, dpSC_insar, dpSC_gps]
% Adjust the indices based on your actual pFromV output structure
pressure_samples = zeros(4, size(posterior2, 2)); 
optParams3 = zeros(4, 1);

for k = 1:size(posterior2, 2)
    temp_m_insar = get_full_m(optimizedM, posterior2(:, k), true, "insar", true); 
    temp_m_gps = get_full_m(optimizedM, posterior2(:, k), true, "gps", true);
    % Calculate pressure
    % Assuming pFromV returns a vector or struct of pressures
    pressure_samples(1,k) = spheroid_pFromV(temp_m_insar(1:8),0.25, 3.08*10^9,'volume');
    pressure_samples(2,k) = spheroid_pFromV(temp_m_gps(1:8),0.25, 3.08*10^9,'volume');
    pressure_samples(3,k) = spheroid_pFromV(temp_m_insar(9:16),0.25, 3.08*10^9,'volume');
    pressure_samples(4,k) = spheroid_pFromV(temp_m_gps(9:16),0.25, 3.08*10^9,'volume');
end
for i = 1:4
    pressure_samples(i,pressure_samples(i,:) < -50e6) = nan;
    optParams3(i) = spheroid_pFromV([optimizedM(1:7), optParams2(i+2)],0.25, 3.08*10^9,'volume');
end
optimizedM(8) = optParams3(2); optimizedM(16) = optParams2(5);
% Define the full list of parameters to plot in order (4x4 grid)
finalParamNames = { ...
  '$\Delta V_{\mathrm{HMM}}^{\mathrm{InSAR}}$', '$\Delta V_{\mathrm{HMM}}^{\mathrm{GPS}}$', '$V_{\mathrm{HMM}}$', '$x_{\mathrm{HMM}}$', ...
  '$y_{\mathrm{HMM}}$', '$d_{\mathrm{HMM}}$', '$\alpha_{\mathrm{HMM}}$', '$x_{\mathrm{SC}}$', ...
  '$y_{\mathrm{SC}}$', '$d_{\mathrm{SC}}$', '$\alpha_{\mathrm{SC}}$', '$\phi_{\mathrm{SC}}$', ...
  '$\psi_{\mathrm{SC}}$', '$\Delta V_{\mathrm{SC}}^{\mathrm{InSAR}}$', '$\Delta V_{\mathrm{SC}}^{\mathrm{GPS}}$', '$V_{\mathrm{SC}}$' ...
};

plotScales = [ ...
    1e-9, 1e-9, 1e-9, 1e-3, ...   % Row 1: dP_HMM, dP_HMM, V_HMM, x_HMM
    1e-3, 1e-3, 1,    1e-3, ...   % Row 2: y_HMM, d_HMM, a_HMM, x_SC
    1e-3, 1e-3, 1,    1,    ...   % Row 3: y_SC,  d_SC,  a_SC,  phi_SC
    1,    1e-9, 1e-9, 1e-9  ...   % Row 4: psi_SC, dP_SC, dP_SC, V_SC
];

% Define Source for each slot: 
% 1=Run1(posterior), 2=Run2(posterior2), 3=Calculated(pressure_samples)
% Mapping indices for lookup from the source arrays
sourceMap = [ ...
    2, 3;  2, 4;  2, 1;  1, 3; ... % Row 1 (Note: indices refer to row in source array)
    1, 4;  1, 5;  1, 6;  1, 7; ... % Row 2
    1, 8;  1, 9;  1, 10; 1, 11; ... % Row 3
    1, 12; 2, 5;  2, 6;  2, 2  ... % Row 4
];

% 1 = paramNames, 2 = paramNames2, 3 = paramNames3 (pressures)


figure(5); clf;
tl = tiledlayout(4,4,'Padding','compact', 'TileSpacing','compact');

for i = 1:16
    ax = nexttile;
    
    % Retrieve data based on source map
    srcType = sourceMap(i, 1);
    srcIdx  = sourceMap(i, 2);
    
    if srcType == 1
        data = posterior(srcIdx, :); % From Run 1 (Geometry)
        priorName = paramNames{srcIdx};
        mle = optParams(srcIdx);
        lb_chosen = lb(srcIdx);
        ub_chosen = ub(srcIdx);
    elseif srcType == 2
        data = posterior2(srcIdx, :); % From Run 2 (Volumes)
        priorName = paramNames2{srcIdx};
        mle = optParams2(srcIdx);
        lb_chosen = lb2(srcIdx);
        ub_chosen = ub2(srcIdx);
    else
        data = pressure_samples(srcIdx, :); % From Calculated Pressures
        priorName = paramNames3{srcIdx};
        mle = optParams3(srcIdx);
        lb_chosen = lb3(srcIdx);
        ub_chosen = ub3(srcIdx);
    end
    
    % Create Histogram
    s = plotScales(i);
    [counts, edges] = histcounts(data(~isnan(data)), 50, 'Normalization', 'pdf');
    binCenters = (edges(1:end-1) + diff(edges)/2) * s;
    counts = counts / s;
    
    % Plot Posterior
    hPost = plot(binCenters, counts, 'Color',[0 0.4470 0.7410], 'LineWidth', 3);
    hold on;
    
    % Plot MLE / Mean (For pressure/vol, use mean of samples; for geom use optParams)
    if srcType == 1
        val_mle = optParams(srcIdx) * s;
    else
        val_mle = mean(data) * s; % Use mean for derived parameters
    end
    xline(mle*s, '--r', 'LineWidth', 2);
    
    % Plot Priors
    % Extend prior beyond posterior range by 10% each side
    span = edges(end) - edges(1);
    xGrid = linspace(lb_chosen, ub_chosen, 400)';
    name = finalParamNames{i};
    p_prior = pdf(paramDists.(priorName).dist, xGrid);
    
    % Rescale prior 
    xGrid = xGrid * s;
    p_prior = p_prior / s;
    
    % Prior PDF in dark gray dashed 
    % if(i ~= 14 && i ~= 15 && i ~= 16)
        hPrior = plot(xGrid, p_prior, '--', 'Color',[0.3 0.3 0.3], 'LineWidth',2.5);
    % end

    xlim([lb_chosen * s, ub_chosen * s]);
    % Formatting
    grid on;
    title(finalParamNames{i}, 'Interpreter','latex', 'FontSize', 16);
    axis square;
    
    % Print CI stats
    ci = prctile(data, [5 95]) * s;
    fprintf('%s: 90%% CI [%.2f, %.2f]\n', finalParamNames{i}, ci(1), ci(2));
end

% Add Legend
if saveFigs; exportgraphics(tl, './PaperFigs/combined_hist.png', 'Resolution', 300); end


%% Plot correlation btwn parameters
% plotParamNames_nounit = {
%   '$\Delta p_{\mathrm{HMM}}^{\mathrm{InSAR}}$', ...
%   '$\Delta p_{\mathrm{HMM}}^{\mathrm{GPS}}$', ...
%   '$V_{\mathrm{HMM}}$', ...
%   '$x_{\mathrm{HMM}}$', ...
%   '$y_{\mathrm{HMM}}$', ...
%   '$d_{\mathrm{HMM}}$', ...
%   '$\alpha_{\mathrm{HMM}}$', ...
%   '$x_{\mathrm{SC}}$', ...
%   '$y_{\mathrm{SC}}$', ...
%   '$d_{\mathrm{SC}}$', ...
%   '$\alpha_{\mathrm{SC}}$', ...
%   '$\phi_{\mathrm{SC}}$', ...
%   '$\psi_{\mathrm{SC}}$', ...
%   '$\Delta p_{\mathrm{SC}}^{\mathrm{InSAR}}$', ...
%   '$\Delta p_{\mathrm{SC}}^{\mathrm{GPS}}$', ...
%   '$V_{\mathrm{SC}}$',
% };
% % assume posterior is already (nParams × nSamples)
% nparams = size(posterior,1);
% f = figure(12); clf;
% tl = tiledlayout(nparams,nparams,'Padding','tight','TileSpacing','tight');
% set(gcf,'Color','w');  % white background
% set(f, 'Units', 'pixels'); 
% set(f, 'Position', [100, 100, 1800, 1600]);
% 
% numBins = 30;
% 
% for i = 1:nparams
%     for j = 1:nparams
%         ax = nexttile(tl);
%         hold(ax,'on');
%         set(ax,'FontSize',8);
%         if i > j
%             % 2D density in lower triangle
%             postx = posterior(j,:);
%             posty = posterior(i,:);
%             [N, edgesX, edgesY] = histcounts2(postx,posty,numBins);
%             centersX = (edgesX(1:end-1)+edgesX(2:end))/2;
%             centersY = (edgesY(1:end-1)+edgesY(2:end))/2;
%             [X,Y] = meshgrid(centersX,centersY);
%             contourf(ax, X, Y, N.', 10, 'LineColor','none');
%             colormap(ax,'hot');
%             grid(ax,'on');
%             ax.YAxis.Visible = 'off';
%         elseif i == j
%             % 1D histogram on diagonal
%             histogram(ax, posterior(i,:), numBins, 'FaceColor',[.2 .6 .5]);
%             % title(ax, plotParamNames_nounit{i}, ...
%             %       'Interpreter','latex','FontSize',16);
%             text(ax, 0.5, 1.45, plotParamNames_nounit{i}, ...
%                 'Units', 'normalized', ... 
%                 'HorizontalAlignment', 'center', ...
%                 'Interpreter', 'latex', ...
%                 'FontSize', 28, ... % You can now make this huge without squishing
%                 'FontWeight', 'bold');
%             ax.YAxis.Visible = 'off';
%         else
%             % blank above diagonal
%             axis(ax,'off');
%             continue
%         end
% 
%         % restore tick labels everywhere
%         ax.XAxis.Visible = 'off';
% 
% 
% 
%         % only bottom row: add xlabel
%         if i == nparams
%             ax.XLabel.String = plotParamNames_nounit{j};
%             ax.XLabel.Interpreter = 'latex';
%             ax.XLabel.FontSize = 20;
%             ax.XAxis.Visible = 'on';
%         else
%             ax.XLabel.String = '';
%         end
% 
%         % only first column: add ylabel
%         if j == 1
%             ax.YLabel.String = plotParamNames_nounit{i};
%             ax.YLabel.Interpreter = 'latex';
%             ax.YLabel.FontSize = 20;
%             ax.YAxis.Visible = 'on';
%         else
%             ax.YLabel.String = '';
%         end
%     end
% end
% 
% sgtitle(tl, "2D Density (lower triangle) & Histograms (diag.) for MCMC", ...
%         'FontWeight','normal', 'FontSize', 30);
% 
% if saveFigs
%     exportgraphics(tl, './PaperFigs/corr_plot.png','Resolution',500);
% end
% 
% clearvars i j options resnorm resdpos dsize dangle A b res;

%% Creating greens functions using the new parameters for M
[gHMM, gSC] = creategreens(optimizedM(1:8), optimizedM(9:end));
[gTiltHMM, gTiltSC] = createtiltgreens(optimizedM(1:8), optimizedM(9:end), 0, false);

%% Setting up offsets to optimize for
nanstatbeginning = [isnan(ux(:, 1)); isnan(uy(:, 1)); isnan(uz(:, 1))];

%% Calculate the pressure at each time step
% Setting up green's functions and formatting displacements/tilts properly
gHMMflat = gHMM';
gHMMflat = gHMMflat(:);
gSCflat = gSC';
gSCflat = gSCflat(:);

ntime = max(size(u(1,1,:)));

tiltx = tilts(:, 1); 
tilty = tilts(:, 2); 

i = 1;
dispstd = [invStdPWRL(1) .* ones(1, length(ux(:, i))), invStdPWRL(2) .* ones(1, length(ux(:, i))), invStdPWRL(3) .* ones(1, length(ux(:, i))), 1, 1];

dp_weight = 1e6;
temp = TimeDependentLSQtilt(gHMMflat, gSCflat, gTiltHMM, gTiltSC, ux, uy, uz, tiltx, tilty, ...
    dispstd, GPSNameList, rw_stddev, dp_weight, true);
dp = real([temp(1:length(tiltx)); temp((length(tiltx) + 1):(2*length(tiltx)))])';
offsets = temp((2*length(tiltx) + 1):end);
ones_temp = zeros(1, length(nanstatbeginning));
ones_temp(nanstatbeginning) = offsets;
offsets = ones_temp;
nanstatbeginning = nanstatbeginning(1:nstat);
clear temp i ones_temp;

%% Modify ux, uy, uz to reflect solved offsets
offsets = reshape(offsets, [], 3);
% ux(nanstatbeginning, :) = ux(nanstatbeginning, :) + offsets(:, 1);
% uy(nanstatbeginning, :) = uy(nanstatbeginning, :) + offsets(:, 2);
% uz(nanstatbeginning, :) = uz(nanstatbeginning, :) + offsets(:, 3);

% Fill outliers, which occur due to gaps in data throwing off LSQ suddenly
% Only fill outliers based on SC data bc HMM has many "outliers" due to
% pressure spikes
[dp(:, 2), TF] = filloutliers(dp(:, 2), "makima", "movmedian", 100, 1);
dp(TF, 1) = NaN;
dp(:, 1) = fillmissing(dp(:, 1), "makima", 1);

%% Error analysis

N_draws = 10;
N_noise = 5;

disp("Getting errors...")
[dp_low, dp_high, u_low, u_high] = GetErrors(N_draws, N_noise, posterior, posterior2, pressure_samples, paramDists, ntime, ux, uy, uz, tiltx, tilty, dispstd, ...
    GPSNameList, rw_stddev, dp_weight, taiyi_parameters, npitloc, invStdPWRL, nanstatbeginning, paramNames, paramNames2, paramNames3, optimizedM);

%% Creating a matrix to store all predicted displacements + tilts in var called usim
%all errors are stored in u_error_low/high
usim = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));
u_error_low = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));
u_error_high = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));
uobs = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));

for i = 1:max(size(dp))
    % **** SETTING THE SIMULATED DISPLACEMENTS ****
    usim(i, :, 1:nstat) =  (gHMM .* dp(i, 1)) + (gSC .* dp(i, 2));
    usim(i, :, nstat + 1) = [gTiltHMM .* dp(i, 1) + gTiltSC .* dp(i, 2), 0];
end


%% Get chi^2 statistic

% Making convenient matrix for comparing observed and predicted
% displacements
uobs(:, :, 1:nstat) = permute(u(1:nstat, :, :), [3, 2, 1]);
uobs(:, :, nstat + 1) = tilts;

residuals = usim - uobs;

chi2 = 0;

for i = 1:size(usim, 2) % Component of GPS station
    for j = 1:nstat % GPS station number (tilt is the 15th element)
        if(j < nstat + 1); stddev = 1/invStdPWRL(i);
        else; stddev = tiltstd; end

        for k = 1:size(usim, 1) % Number of timesteps for each station
            % Compute contribution to chi2
            s = ((uobs(k, i, j) - usim(k, i, j))^2)/(stddev^2);
            if(~isnan(s))
                chi2 = chi2 + ((uobs(k, i, j) - usim(k, i, j))^2)/(stddev^2);
            end
        end
    end
end

DOF = nnz(~isnan(uobs)) - length(t)*2; % Number of non-nan observations - free params

disp("Reduced chi^2 (no tilt) = " + chi2/DOF);

%% Extract collapse amplitudes
collapset = D.gps_info.t_events(26:end);
collapset = decyear(datetime(collapset, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy'));
% [ampHMM, ampSC] = ExtractCollapseAmplitude([dp(:, 1)'; dp(:, 2)'], t, collapset, t(3) - t(2));


%% Plot shear stress change
nu = 0.25;
mu = 3.08*10^9;
R = 900; % m
H = 800; % m
g = 9.8; % m/s^2
rho_c = 2500; % kg/m^3

% Using equation 1 (assuming dv/dt = 0) in PNAS paper:
% dp must be rescaled by the greens function pressure drop but that's
% already dealt with for the error bnds.
tau_pnas = - R/(2 * H) * dp(:,1) * optimizedM(8); 
% Must convert the bnds from MPa to Pa
tau_pnas_low = -R/(2 * H) * dp_low(:,1)*1e6;
tau_pnas_high = -R/(2 * H) * dp_high(:,1)*1e6;
%
% figure(3); clf;
% hold on;
% plot(tau_pnas .* 1e-6);
% legend();
% ylabel("Average tau (MPa)")
% xlabel("Time"); axis square;
%%
makeplots(x, y, GPS_llh, u, u1d, ux, uy, uz, u_low, u_high, tiltx, tilty, usim, t, nanstatend, ...
    nanstatbeginning, finalindex, collapset, dp, dp_low, dp_high, tau_pnas, tau_pnas_low, tau_pnas_high, optParams, optimizedM, ...
    GPSNameList, gTiltHMM, gTiltSC, xtilt, ytilt, tiltreduced, radscale, coast_new, taiyi_parameters, 3, ntrials, offsets, saveFigs);

%% Insar plotting
% Set if we want to plot the ascending or descending track with insarmode
insarmode = "asc";
% make predicted insar data
if(insarmode == "asc")
    ind = (insar_lengths(1) + 1):sum(insar_lengths);
    insaru_pred = real(spheroid(optimizedM(1:8), [insarx(ind); insary(ind); zeros(size(insarx(ind)))], 0.25, 3.08*10^9) + ...
        spheroid(optimizedM(9:end), [insarx(ind); insary(ind); zeros(size(insarx(ind)))], 0.25, 3.08*10^9));
    insaru_pred = insaru_pred' * look(:,2);
else
    ind = 1:(insar_lengths(1));
    insaru_pred = real(spheroid(optimizedM(1:8), [insarx(ind); insary(ind); zeros(size(insarx(ind)))], 0.25, 3.08*10^9) + ...
        spheroid(optimizedM(9:end), [insarx(ind); insary(ind); zeros(size(insarx(ind)))], 0.25, 3.08*10^9));
    insaru_pred = insaru_pred' * look(:,1);
end

cLimits = [-1.55, 0.55];
opacity = 0.7;
cmap = turbo;
% Set if x and y axis labels are on
xon = true; yon = false;
% Set if colorbar is on
con = true;

plot_insar_new(insarx(ind), insary(ind), insaru_pred', block_size(ind), look, x, y, u1d, u1d, xon, yon, con, ...
    31, GPSNameList, optimizedM, coast_new,cLimits, opacity, saveFigs, insarmode);
%  - insaru_pred' insaru_full(ind)