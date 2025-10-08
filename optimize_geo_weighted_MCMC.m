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
% taiyi_parameters = [1600.79, 914.47, 90, 0, 0.46e3 + npitloc(1), 0.35e3 + npitloc(2), -2.18e3, -4e7, ... 
%     277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, -10e6];
% Use taiyi + kyle for HMM 
taiyi_parameters = [1600.79, 914.47, 90, 0, 50, 200, -2.18e3, -4e7, ... 
     277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, -10e6];

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
% ["dpHMM_insar", "dpHMM_gps", "volHMM", 'xHMM', 'yHMM', "zHMM" (from top of spheroid), 
% xSC, ySC, dSC, "SC aspect ratio", "dip", strike, "dpSC_insar", "dpSC_gps", volSC]
% yHMM lb used to be -0.5e3 + npitloc(2)
% Old:
lb = [-5e7, -5e7, 1e8, -5e2, -2e3 + npitloc(2), -2e3, 0.8, ...
    0.7e3, -2.9e3, -4.6e3, 0.05, 40, 100, -5e7, -5e7, 2.4e9]; 
ub = [1e6, 1e6, 3e10, 0.5e3 + npitloc(1), 0.5e3 + npitloc(2), -3e2, 2.0, ...
     1.8e3, -2.0e3, -2.9e3, 1, 90, 180, 1e6, 1e6, 20e9]; 

% Inferred from kyle / taiyi:
% lb = [-5e7, -5e7, 1e8, -100, 0, -1.3e3, ...
%     -2.5e3 + npitloc(1), -3.4e3 + npitloc(2), -4.7e3, 0.05, 40, 100, -5e7, -5e7, 3e9];
% 
% ub = [1e6, 1e6, 1e10, 200, 400, -200, ...
%     2.5e3 + npitloc(2), -1e3 + npitloc(2), -2.7e3, 0.3, 90, 180, 1e6, 1e6, 20e9];

% Save figures for export
saveFigs = false;

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

%% MCMC Static inversion
paramNames = {'dpHMM_insar', 'dpHMM_gps', 'volHMM', 'xHMM', 'yHMM', 'dHMM', 'alphaHMM' ...
    'xSC', 'ySC', 'dSC', 'alphaSC', 'dipSC', 'strikeSC', 'dpSC_insar', 'dpSC_gps', 'volSC'};
ntrials = 1e5; % Customize to get convergence

% Testing GPS and prior weights.
gps_weights = linspace(4e1, 8e1, 10);
prior_weights = linspace(3e2, 1e3, 10);
gps_weight = 6.7e1; % Optimal weight based on L curve
prior_weight = 5.3e2; % Optimal weight based on L curve

delete(gcp('nocreate'));

% Set up log likelihood lists for L curve creation. 
gps_l2s = zeros(1,length(gps_weights));
insar_l2s = zeros(1,length(gps_weights));
posteriors_list = zeros(length(paramNames), ntrials - 4e3 + 1, length(gps_weights));
optParams_list = zeros(length(paramNames), length(gps_weights));
l_curve_points = zeros(3,length(gps_weights));

runMCMC = true;

% For loop to create L curve
% for i = 1:length(gps_weights)

prior_params = taiyi_parameters;
prior_params(7) = prior_params(7) + prior_params(1);
if(runMCMC)
    [optParams, posterior, gps_l2, insar_l2, prior_l2] = optimize_SC_MCMC(prior_params, lb, ub, xopt, ...
        yopt, zopt, u1d', insarx, insary, insaru_full', look, insar_lengths, sparse(cinv_full), daily_inv_std, ...
        nanstatend, ntrials, gps_weight, prior_weight, paramNames, saveFigs);
    optParams = real(optParams');
    save Data/MCMC_1e6_SCvol.mat optParams posterior gps_l2 insar_l2 prior_l2;
else
    load Data/MCMC_vars_1e6_allparams_nodpineq.mat;
end

% More L curve stuff
%     l_curve_points(1,i) = gps_l2;
%     l_curve_points(2,i) = insar_l2;
%     l_curve_points(3,i) = prior_l2;
% 
%     posteriors_list(:, :, i) = posterior;
%     optParams_list(:, i) = optParams;
% end


% Save / load L curve data
% save Data/l_curve_data_prior.mat l_curve_points;
% load Data/l_curve_data_prior.mat;
% posterior = squeeze(posteriors_list(:, :, i));
% optParams = optParams_list(:, i)';
% disp(gps_weights(end));


% Get the full geometry parameters based on the optimization results:
disp("GPS L2: " + gps_l2 + " InSAR L2: " + insar_l2);
optimizedM = get_full_m(taiyi_parameters, optParams, true, "insar");

% end

%%
% Show parameter table
disp(array2table(optParams(:).', 'VariableNames',paramNames));
% Make sure parameters are real
optParams = real(optParams);

% % L curve plotting 
% % Prepare L-curve data, GPS vs. insar
% % xl = l_curve_points(1,:);
% % yl = l_curve_points(2,:);
% % xl = zeros(100, 1);
% % yl = zeros(100, 1);
% 
% % L-curve data, GPS + insar vs. prior
% xl = l_curve_points(1,:) + l_curve_points(2,:);
% yl = abs(l_curve_points(3,:));
% 
% % Scatter, mapping color to gps_weights
% l_curve = figure(7); clf;
% scatter_handle = scatter(xl, yl, 400, prior_weights, 'filled');
% 
% % Choose colormap and add colorbar
% colormap(parula);
% c = colorbar;               
% c.Label.String = 'Prior Weight';
% c.FontSize = 24;
% 
% % (Optional) fix the color limits if you want to highwlight a subset:
% % caxis([min(gps_weights), max(gps_weights)]);
% 
% % Labels & title
% % xlabel("GPS L2");
% % ylabel("InSAR L2");
% xlabel("GPS L2 Norm + InSAR L2 Norm", 'FontSize', 24);
% ylabel("Prior Log Likelihood", 'FontSize', 24);
% title("L-curve colored by prior weight", 'FontSize', 30);
% axis square;
% ax = gca;
% ax.FontSize = 24;
% 
% % Set up interactive data cursor
% dcm = datacursormode(gcf);
% set(dcm, 'UpdateFcn', @(~, event_obj) ...
%     sprintf('GPS L2: %.2e\nInSAR L2: %.2e\nGPS Weight: %.1e', ...
%             event_obj.Position(1), ...
%             event_obj.Position(2), ...
%             gps_weights(event_obj.DataIndex)));
% 
% exportgraphics(l_curve, './PaperFigs/l_curve_prior.png', 'Resolution', 500);


%% Plot histogram of each parameter
plotParamNames = {
  '$\Delta p_{\mathrm{HMM}}^{\mathrm{InSAR}}\ (\mathrm{MPa})$', ...
  '$\Delta p_{\mathrm{HMM}}^{\mathrm{GPS}}\ (\mathrm{MPa})$', ...
  '$V_{\mathrm{HMM}}\ (\mathrm{km}^3)$', ...
  '$x_{\mathrm{HMM}}\ (\mathrm{km})$', ...
  '$y_{\mathrm{HMM}}\ (\mathrm{km})$', ...
  '$d_{\mathrm{HMM}}\ (\mathrm{km})$', ...
  '$\alpha_{\mathrm{HMM}}$', ...
  '$x_{\mathrm{SC}}\ (\mathrm{km})$', ...
  '$y_{\mathrm{SC}}\ (\mathrm{km})$', ...
  '$d_{\mathrm{SC}}\ (\mathrm{km})$', ...
  '$\alpha_{\mathrm{SC}}$', ...
  '$\phi_{\mathrm{SC}}\ (^\circ)$', ...
  '$\psi_{\mathrm{SC}}\ (^\circ)$', ...
  '$\Delta p_{\mathrm{SC}}^{\mathrm{InSAR}}\ (\mathrm{MPa})$', ...
  '$\Delta p_{\mathrm{SC}}^{\mathrm{GPS}}\ (\mathrm{MPa})$', ...
  '$V_{\mathrm{SC}}\ (\mathrm{km}^3)$',
};

% Scale units to appropriate factor
unitScaling = [1e-6, 1e-6, 1e-9, 1e-3, 1e-3, 1e-3, 1, ...
    1e-3, 1e-3, 1e-3, 1, 1, 1, 1e-6, 1e-6, 1e-9];

figure(5);
clf;
tl = tiledlayout(4,4,'Padding','compact', 'TileSpacing','compact');
priormeans = get_full_m(taiyi_parameters, [], false, "insar");
load Data/paramDists.mat;

for i = 1:length(paramNames)
    ax = nexttile;
    numBins = 100;
    [counts, edges] = histcounts(posterior(i, :), numBins, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;

    % Rescale for unit conversion
    s = unitScaling(i);
    binCenters = binCenters * s;
    counts = counts / s;
    
    % Posterior PDF line in blue 
    hPost = plot(binCenters(2:end), counts(2:end), 'Color',[0 0.4470 0.7410], 'LineWidth',4);
    hold on;
    
    % MLE estimate as red dashed line 
    hMLE = xline(optParams(i) * s, '--r', 'LineWidth',2.5);
    
    % Extend prior beyond posterior range by 10% each side
    span = edges(end) - edges(1);
    xGrid = linspace(lb(i), ub(i), 400)';
    name = plotParamNames{i};
    p_prior = pdf(paramDists.(paramNames{i}).dist, xGrid);
    
    % Rescale prior 
    xGrid = xGrid * s;
    p_prior = p_prior / s;
    
    % Prior PDF in dark gray dashed 
    if(i ~= 13 && i ~= 14)
        hPrior = plot(xGrid, p_prior, '--', 'Color',[0.3 0.3 0.3], 'LineWidth',2.5);
    end
    
    xlim([lb(i)*s, ub(i)*s]);
    % Aesthetics 
    grid on;
    title(name, 'Interpreter','latex', 'FontSize', 30, "FontWeight","bold");
    set(gca, 'YTickLabel', [], 'FontSize', 20);
    
    % compute and print the 90% prior interval
    cdf_values = cumtrapz(xGrid, p_prior);
    [cdf_u, ia] = unique(cdf_values);
    x_u = xGrid(ia);
    lower_bound = interp1(cdf_u, x_u, 0.1);
    upper_bound = interp1(cdf_u, x_u, 0.9);
    fprintf('%s 90%% CI: [%.1e, %.1e], MLE = %.2e\n', ...
            name, lower_bound, upper_bound, optParams(i)*s);
    axis square
end

% Manually add legend to empty tile (tile 16 in a 4x4 grid)
% axLegend = nexttile(15); % explicitly go to the 16th tile
% axis(axLegend, 'off');   % make this tile invisible

% Create legend manually
% leg = legend(axLegend, [hPost, hMLE, hPrior], ...
%        {'Posterior PDF', 'MLE estimate', 'Prior PDF'}, ...
%        'FontSize', 24, 'Box', 'on');
% 
% set(leg, 'Units', 'normalized', 'Position', [0.793, 0.212, 0.172, 0.1]);
% hold off;

if saveFigs
    exportgraphics(tl, './PaperFigs/geo_hist.png', 'Resolution', 500);
end


gps_l2s(i) = gps_l2;
insar_l2s(i) = insar_l2;

disp("MCMC done");

save("./Figures/optimizedM.mat");

%% Plot correlation btwn parameters
% assume posterior is already (nParams × nSamples)
nparams = size(posterior,1);
figure(12); clf;
tl = tiledlayout(nparams,nparams,'Padding','tight','TileSpacing','tight');
set(gcf,'Color','w');  % white background

numBins = 30;

for i = 1:nparams
    for j = 1:nparams
        ax = nexttile(tl);
        hold(ax,'on');
        set(ax,'FontSize',8);
        if i > j
            % 2D density in lower triangle
            postx = posterior(j,:);
            posty = posterior(i,:);
            [N, edgesX, edgesY] = histcounts2(postx,posty,numBins);
            centersX = (edgesX(1:end-1)+edgesX(2:end))/2;
            centersY = (edgesY(1:end-1)+edgesY(2:end))/2;
            [X,Y] = meshgrid(centersX,centersY);
            contourf(ax, X, Y, N.', 10, 'LineColor','none');
            colormap(ax,'hot');
            grid(ax,'on');
            ax.YAxis.Visible = 'off';
        elseif i == j
            % 1D histogram on diagonal
            histogram(ax, posterior(i,:), numBins, 'FaceColor',[.2 .6 .5]);
            title(ax, plotParamNames{i}, ...
                  'Interpreter','latex','FontSize',16);
            ax.YAxis.Visible = 'off';
        else
            % blank above diagonal
            axis(ax,'off');
            continue
        end

        % restore tick labels everywhere
        ax.XAxis.Visible = 'off';
        
        

        % only bottom row: add xlabel
        if i == nparams
            ax.XLabel.String = plotParamNames{j};
            ax.XLabel.Interpreter = 'latex';
            ax.XLabel.FontSize = 12;
            ax.XAxis.Visible = 'on';
        else
            ax.XLabel.String = '';
        end

        % only first column: add ylabel
        if j == 1
            ax.YLabel.String = plotParamNames{i};
            ax.YLabel.Interpreter = 'latex';
            ax.YLabel.FontSize = 12;
            ax.YAxis.Visible = 'on';
        else
            ax.YLabel.String = '';
        end
    end
end

sgtitle(tl, "2D Density (lower triangle) & Histograms (diag.) for MCMC", ...
        'FontWeight','normal', 'FontSize', 30);

if saveFigs
    exportgraphics(tl, './PaperFigs/corr_plot.png','Resolution',500);
end

clearvars i j options resnorm resdpos dsize dangle A b res;

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

N_draws = 20;
N_noise = 5;

disp("Getting errors...")
[dp_low, dp_high, u_low, u_high] = GetErrors(N_draws, N_noise, posterior, paramDists, ntime, ux, uy, uz, tiltx, tilty, dispstd, ...
    GPSNameList, rw_stddev, dp_weight, taiyi_parameters, npitloc, invStdPWRL, nanstatbeginning, paramNames, optParams);

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

%% Make synthetic insar data based on optimized geometry
ind = 1:(insar_lengths(1));
% ind = (insar_lengths(1) + 1):sum(insar_lengths);
insaru_pred = spheroid(optimizedM(1:8), [insarx(ind); insary(ind); zeros(size(insarx(ind)))], 0.25, 3.08*10^9) + ...
    spheroid(optimizedM(9:end), [insarx(ind); insary(ind); zeros(size(insarx(ind)))], 0.25, 3.08*10^9);
insaru_pred = insaru_pred' * look(:,1);

%% Plot shear stress change
nu = 0.25;
mu = 3.08*10^9;
R = 700; % m
H = 500; % m
g = 9.8; % m/s^2
rho_c = 2500; % kg/m^3

% Using equation 1 (assuming dv/dt = 0) in PNAS paper:

tau_pnas = R * rho_c * g / 2 - R/(2 * H) * dp(:,1) * optimizedM(8);
%
figure(3); clf;
hold on;
plot(tau_pnas .* 1e-6);
legend();
ylabel("Average tau (MPa)")
xlabel("Time"); axis square;
%%
makeplots(x, y, GPS_llh, u, u1d, ux, uy, uz, u_low, u_high, insarx(ind), insary(ind), insaru_full(ind), insaru_pred, block_size(ind), look(:,1), tiltx, tilty, usim, t, nanstatend, ...
    nanstatbeginning, finalindex, collapset, dp, dp_low, dp_high, tau_pnas, optParams, optimizedM, GPSNameList, gTiltHMM, ...
    gTiltSC, xtilt, ytilt, tiltreduced, radscale, coast_new, taiyi_parameters, 3, ntrials, offsets, saveFigs);

%% Insar plotting
cLimits = [-1.55, 0.55];
opacity = 0.7;
cmap = turbo;
plot_insar_new(insarx(ind), insary(ind), insaru_full(ind) - insaru_pred', block_size(ind), look(:,1), x, y, u1d, u1d, nanstatend, ...
radscale, 31, GPSNameList, optimizedM, coast_new,cLimits, opacity, saveFigs);
%  - insaru_pred' insaru_full(ind)