function [dp_low, dp_high, u_high, u_low] = GetErrors(N_draws, N_noise, posterior, paramDists, ntime, ...
    ux, uy, uz, tiltx, tilty, ...
    dispstd, GPSNameList, rw_stddev, dp_weight, taiyi_parameters, npitloc, gps_sigma, nanstatbeginning, paramNames, optParams)
%% First construct the covariance matrices to generate added nosie
%%% Delete dp, optimizedM when it works
N = length(tiltx);          % number of measurements
C_rw = zeros(N,N);
parfor i = 1:N
    for j = 1:N
        C_rw(i,j) = rw_stddev^2 * min(i,j);
    end
end
gps_sigma = 1./gps_sigma;

%% Use posterior distribution to generate a list of many candidate geometries
dtheta = 0;
full_param_names = {'dpHMM_insar', 'dpHMM_gps','volHMM', 'alphaHMM', 'xHMM', 'yHMM', 'dHMM', ...
    'alphaSC', 'dipSC', 'strikeSC', 'xSC', 'ySC', 'dSC', 'dpSC_insar', 'dpSC_gps', 'mu'}; 
full_posterior = zeros(length(full_param_names), size(posterior,2));

% Excluding dip HMM, strike HMM, SC volume % ADD MU

% populate the posterior distribution with each relevant metric
for i = 1:size(full_posterior, 1)
    idx = find(strcmp(full_param_names{i}, paramNames));
    if ~isempty(idx)
        full_posterior(i, :) = posterior(idx, :);
    else
        % xlims = paramDists.(full_param_names{i}).xlim;
        % xsamples = linspace(xlims(1), xlims(2), length(posterior));
        % pdf_vals = pdf(paramDists.(full_param_names{i}).dist, xsamples);
        full_posterior(i,:) = random(paramDists.(full_param_names{i}).dist, ...
                             1, size(posterior,2));
    end
    
end
% Now convert full_posterior to a posterior with all the params for
% spheroid.m
complete_posterior = zeros(16, size(posterior,2));
% Vert and horiz semi diam
HMM_vol = full_posterior(3,:);
alpha_HMM = full_posterior(4,:);
opt_vert_sd = (3/(4*pi) .* HMM_vol .* (alpha_HMM.^2)).^(1/3);
opt_horiz_sd = opt_vert_sd./(alpha_HMM);
optimizedM = get_full_m(taiyi_parameters, optParams, true, "insar");

% Check one parameterameter sampling correctly
complete_posterior(1:8,:) = optimizedM(1:8)' * ones(1, length(posterior));

complete_posterior(1,:) = opt_vert_sd;
complete_posterior(2,:) = opt_horiz_sd;
complete_posterior(3,:) = 90 .* ones(1, length(posterior));
complete_posterior(4,:) = zeros(1, length(posterior));
complete_posterior(5,:) = full_posterior(5,:);
complete_posterior(6,:) = full_posterior(6,:);

% Must make the z coordinate relative to the top of the sphere - surface
complete_posterior(7,:) = full_posterior(7,:); % -full_posterior(5,:) + taiyi_parameters(1); now this is already relative to top of sphere
complete_posterior(8,:) = full_posterior(1, :); % dpHMM take insar
mu_dist = full_posterior(end,:); % NOTE: mu is in units log10(Pa)

SC_vol = 2.5e9;
% REVERT THIS BACK
alpha_SC = full_posterior(8,:); 
% optParams(3) * ones(1, length(full_posterior(7,:))); 
opt_vert_sd = (3/(4*pi) .* SC_vol .* (alpha_SC.^2)).^(1/3);
opt_horiz_sd = opt_vert_sd./(alpha_SC);

% Check that one parameter is sampling correctly
% complete_posterior(9:16,:) = optimizedM(9:16)' * ones(1,length(posterior));
complete_posterior(9,:) = opt_vert_sd;
complete_posterior(10,:) = opt_horiz_sd;
complete_posterior(11,:) = full_posterior(9, :);
complete_posterior(12,:) = full_posterior(10, :);
complete_posterior(13,:) = full_posterior(11,:);
complete_posterior(14,:) = full_posterior(12,:); 
complete_posterior(15,:) = full_posterior(13,:);
complete_posterior(16,:) = full_posterior(14, :);


figure(14);
for p = 1:length(full_param_names)
    subplot(4, 4, p);
    histogram(complete_posterior(p, :));
    xline(optimizedM(p), 'color', 'red');
    title("Param number = " + p);
end

% randIdx = randi(size(complete_posterior, 2), [1, N_draws]);

geo_samples = zeros(size(complete_posterior, 1), N_draws);
mu_samples = zeros(1, N_draws);
i = 1;
while i < N_draws + 1
    randIdx = randi(size(complete_posterior, 2), [1, 1]);
    % Check that depth is a reasonable value
    % if(complete_posterior(7, randIdx) + abs(complete_posterior(1, randIdx)) < taiyi_parameters(7) + abs(taiyi_parameters(1)))
    geo_samples(:,i) = complete_posterior(:, randIdx);
    % Add offset to top of HMM to bring center back to center of volume
    geo_samples(7,i) = geo_samples(7,i) - geo_samples(1);
    mu_samples(i) = mu_dist(randIdx);
    i = i + 1;
    % end
end

% Include optimal geometry in one of the samples
geo_samples(:, 1) = optimizedM;
dp_dist = zeros(N_draws, N_noise, ntime*2);
ux_dist = zeros(N_draws, N_noise, size(ux, 1), size(ux, 2));
uy_dist = zeros(N_draws, N_noise, size(ux, 1), size(ux, 2));
uz_dist = zeros(N_draws, N_noise, size(ux, 1), size(ux, 2));
tiltx_dist = zeros(N_draws, N_noise, size(tiltx, 1));
tilty_dist = zeros(N_draws, N_noise, size(tiltx, 1));
% h = waitbar(0,'Analyzing errors...');

if isempty(gcp('nocreate'))
    parpool('local');  % or a specific number of workers: parpool('local',4)
end

parfor i = 1:N_draws
    % Format new geometry into optimizedM array
    m_samp = geo_samples(:,i); % get_full_m(taiyi_parameters, geo_samples(:,i), true);
    mu_samp = 10^mu_samples(i);
    % Create new green's functions
    [gHMM_samp, gSC_samp] = creategreens(m_samp(1:8), m_samp(9:end), mu_samp);
    [gTiltHMM_samp, gTiltSC_samp] = createtiltgreens(m_samp(1:8), m_samp(9:end), dtheta, false, mu_samp);
    gHMMflat_samp = gHMM_samp';
    gHMMflat_samp = real(gHMMflat_samp(:));
    gSCflat_samp = gSC_samp';
    gSCflat_samp = real(gSCflat_samp(:));

    temp = TimeDependentLSQtilt(gHMMflat_samp, gSCflat_samp, gTiltHMM_samp, gTiltSC_samp, ux, uy, uz, tiltx, tilty, ...
        dispstd, GPSNameList, rw_stddev, dp_weight, true);

    dp_ideal = temp(1:2*length(tiltx)); %[dp(:, 1), dp(:, 2)];
    offsets = temp((2*length(tiltx) + 1):end);
    offsets = reshape(offsets, [], 3);
    idx_HMM = (1:length(tiltx));
    idx_SC = (length(tiltx) + 1):2*length(tiltx);
    
    %% Now take the inverted pressure history and generate synthetic data
    % Use that synthetic data and add noise. Then re-invert to get the
    % variance of those individual estimates
    ux_pred = gHMM_samp(1,:)' * dp_ideal(1:length(tiltx)) + gSC_samp(1,:)' * dp_ideal(length(tiltx)+1:2*length(tiltx));
    ux_pred(nanstatbeginning, :) = ux_pred(nanstatbeginning, :); % + offsets(:, 1);

    uy_pred = gHMM_samp(2,:)' * dp_ideal(1:length(tiltx)) + gSC_samp(2,:)' * dp_ideal(length(tiltx)+1:2*length(tiltx));
    uy_pred(nanstatbeginning, :) = uy_pred(nanstatbeginning, :); %+ offsets(:, 2);

    uz_pred = gHMM_samp(3,:)' * dp_ideal(1:length(tiltx)) + gSC_samp(3,:)' * dp_ideal(length(tiltx)+1:2*length(tiltx));
    uz_pred(nanstatbeginning, :) = uz_pred(nanstatbeginning, :); %+ offsets(:, 3);

    tilte_pred = gTiltHMM_samp(1) .* dp_ideal(1:length(tiltx)) + gTiltSC_samp(1) .* dp_ideal(length(tiltx)+1:2*length(tiltx));
    tiltn_pred = gTiltHMM_samp(2) .* dp_ideal(1:length(tiltx)) + gTiltSC_samp(2) .* dp_ideal(length(tiltx)+1:2*length(tiltx));
    
    for j = 1:N_noise
        % Generate noise for each data
        rw_east_noise = mvnrnd(zeros(length(C_rw), 1), C_rw);
        rw_north_noise = mvnrnd(zeros(length(C_rw), 1), C_rw);
        gps_e_noise = normrnd(0,gps_sigma(1),[1,length(ux)]);
        gps_n_noise = normrnd(0,gps_sigma(2),[1,length(ux)]);
        gps_u_noise = normrnd(0,gps_sigma(3),[1,length(ux)]);
    
        % Add generated noise to each data
        tilte_noised = tilte_pred + rw_east_noise;
        tiltn_noised = tiltn_pred + rw_north_noise;
        ux_noised = ux_pred + gps_e_noise;
        uy_noised = uy_pred + gps_n_noise;
        uz_noised = uz_pred + gps_u_noise;
        
        ux_dist(i, j, :, :) = ux_noised;
        uy_dist(i, j, :, :) = uy_noised;
        uz_dist(i, j, :, :) = uz_noised;

        tiltx_dist(i, j, :) = tilte_noised;
        tilty_dist(i, j, :) = tiltn_noised;

        % Make NaN values for ux, uy, uz to match original data
        nanx = isnan(ux); nanx(:, 1) = 0;
        nany = isnan(ux); nany(:, 1) = 0;
        nanz = isnan(ux); nanz(:, 1) = 0;

        ux_noised(nanx) = nan; 
        uy_noised(nany) = nan;
        uz_noised(nanz) = nan;
    
        %% Re-invert for pressure of now re-noised synthetic data
        temp = TimeDependentLSQtilt(gHMMflat_samp, gSCflat_samp, gTiltHMM_samp, gTiltSC_samp, ux_noised, uy_noised, uz_noised, ...
             tilte_noised, tiltn_noised, dispstd, GPSNameList, rw_stddev, dp_weight, true);
        
        % Fill outliers
        [temp(idx_SC), TF] = filloutliers(temp(idx_SC), "makima", "movmedian", 100, 1);
        temp(TF) = nan; % Now fill in HMM outliers
        temp(idx_HMM) = fillmissing(temp(idx_HMM), "makima", 1);

        dp_dist(i, j, :) = [temp(idx_HMM) .* m_samp(8)/(1e6), temp(idx_SC) .* m_samp(16)/(1e6)];
    end
end
dp_dist = reshape(dp_dist, N_draws * N_noise, ntime*2);

clear gHMM_samp gSC_samp gTiltHMM_samp gTiltSC_samp temp

%% Now pick the 10th and 90th percentile points for each sample
% Get conf. interval for pressure
pLow  = prctile(dp_dist, 10, 1)';  % 10th percentile across columns
pHigh = prctile(dp_dist, 90, 1)';  % 90th percentile

dp_low = real([pLow(1:length(tiltx)), pLow((length(tiltx) + 1):(2*length(tiltx)))]);
dp_high = real([pHigh(1:length(tiltx)), pHigh((length(tiltx) + 1):(2*length(tiltx)))]);


% Get conf interval for displacements
% assume ux_dist, uy_dist, uz_dist are each [G × H x S × T]
[G, H, S, T] = size(ux_dist);

% 1) Flatten geom+noise dims into one
ux_flat = reshape(ux_dist, H*G, S, T);
uy_flat = reshape(uy_dist, H*G, S, T);
uz_flat = reshape(uz_dist, H*G, S, T);

[G, H, St] = size(tiltx_dist);

% 1) Flatten geom+noise dims into one
tiltx_flat = reshape(tiltx_dist, H*G, St);
tilty_flat = reshape(tilty_dist, H*G, St);

% 2) Compute percentiles along the first dimension of each flat array
%    prctile(...,p,1) returns [1 × S × T], so squeeze to [S × T]
ux_p10 = squeeze( prctile(ux_flat, 10, 1) );  % [S × T]
ux_p90 = squeeze( prctile(ux_flat, 90, 1) );

uy_p10 = squeeze( prctile(uy_flat, 10, 1) );
uy_p90 = squeeze( prctile(uy_flat, 90, 1) );

uz_p10 = squeeze( prctile(uz_flat, 10, 1) );
uz_p90 = squeeze( prctile(uz_flat, 90, 1) );

tiltx_p10 = squeeze( prctile(tiltx_flat, 10, 1) );
tiltx_p90 = squeeze( prctile(tiltx_flat, 90, 1) );
tilty_p10 = squeeze( prctile(tilty_flat, 10, 1) );
tilty_p90 = squeeze( prctile(tilty_flat, 90, 1) );

% 3) Pack into [Stations × Components × Time] form
%    component‐order: 1=East, 2=North, 3=Up
u_low = zeros(S+1, 3, T);
u_high = zeros(S+1, 3, T);

u_low(1:S,1,:) = ux_p10;
u_low(1:S,2,:) = uy_p10;
u_low(1:S,3,:) = uz_p10;

u_low(end,1,:) = tiltx_p10;
u_low(end,2,:) = tilty_p10;

u_high(1:S,1,:) = ux_p90;
u_high(1:S,2,:) = uy_p90;
u_high(1:S,3,:) = uz_p90;
u_high(end,1,:) = tiltx_p90;
u_high(end,2,:) = tilty_p90;


%% Make a plot of all the diff pressure histories
figure(1);
clf;
% inds = 1:length(tiltx);%(length(tiltx)+1):2*length(tiltx);
% for i = 2:(N_draws*N_noise)
%     plot(dp_dist(i, inds), 'HandleVisibility','off');
%     hold on;
% end

% plot(dp_dist(1, inds), "LineWidth", 4, "DisplayName", "MLE");
% plot(dp_low(:, 1), "LineWidth", 6, "DisplayName", "10th percentile");
% plot(dp_high(:, 1), "LineWidth", 6, "DisplayName", "90th percentile");
plot(tilty);
hold on;
plot(squeeze(u_low(end, 2, :)), "DisplayName", "10th percentile");
plot(squeeze(u_high(end, 2, :)), "DisplayName", "90th percentile");
legend();
end
