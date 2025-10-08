function [optParams, posterior] = optimize_SC_MCMC_par(m_known, lb, ub, xopt, xtilt, yopt, ytilt, zopt, u1d, ...
    invStdPWRL, tiltstd, tiltreduced, nanstatbeginning, ntrials, n_walkers, saveFigs)
% This version uses the parallel gwmcmc sampler instead of the serial mcmc.

    %% Preliminaries
    GPS_std = 1./invStdPWRL;
    vol = (4/3 * pi * m_known(1) * m_known(2)^2);
    priormeans = [vol, m_known(8), m_known(9), m_known(10), m_known(11), m_known(end)];
    paramNames = ["HMM volume", "dpHMM", "vert semi-diameter", "horiz semi-diameter", "dip", "dpSC"];
    
    % Arrange bounds as a matrix with one row per parameter
    rescale_factor = 1./abs(priormeans)';
    bnds = [lb.*rescale_factor'; ub.*rescale_factor']';
    
    % Set up sigma for likelihood weighting.
    sigma = ones(size(u1d));
    for i = 1:3
        sigma(:, i) = sigma(:, 1) * GPS_std(i);
    end
    
    %% Set Up the Ensemble for gwmcmc
    num_params = length(priormeans);  % should be 6
    num_walkers = n_walkers;      % e.g., 12 walkers
    % Create an ensemble by perturbing the initial guess
    minit = repmat(priormeans(:).*rescale_factor, 1, num_walkers) + 0.05 * randn(num_params, num_walkers);
    
    %% Define the Log-Probability Function
    % This anonymous function wraps the model prediction (create_MCMC_data)
    % and computes the log likelihood exactly as in your mcmc routine.
    log_prob = @(m) local_log_prob(m, u1d, priormeans, ntrials, m_known, xopt, xtilt, yopt, ytilt, ...
                                  zopt, invStdPWRL, tiltstd, tiltreduced, nanstatbeginning, sigma, bnds, rescale_factor);
    
    %% Run the Parallel MCMC Sampler
    % We set 'OutputData' to 0 so that the log_prob function only returns one output.
    [models, logP, ~, RunTime, acceptrate] = gwmcmc(minit, {log_prob}, ...
        ntrials, 'Parallel', true, 'NCores', 8, 'OutputData', 1, 'StepSize', 0.92, 'ThinChain', 1);

    
    fprintf('Acceptance ratio: %.2f%%\n', acceptrate * 100);
    
    %% Post-Processing
    % Discard burn-in samples (adjust burn as needed)
    burn = 15;
    models(:,:,1:burn) = [];
    
    % Collapse the ensemble: models is [num_params x num_walkers x num_samples]
    samples = reshape(models, num_params, [])' ./ rescale_factor';
    
    % For optimal parameters we take the sample with the maximum log likelihood.
    L_all = reshape(logP(:, :, (burn+1):end), 1, []);
    [~, idx] = max(L_all);
    optParams = samples(idx, :)';
    posterior = samples;  % all samples after burn-in
    
    %% Plotting (Histograms, evolution, etc.)
    % Histogram plot
    figure(5);
    for i = 1:length(paramNames)
        subplot(320 + i);
        numBins = 100;
        [counts, edges] = histcounts(samples(:,i), numBins, 'Normalization', 'pdf');
        binCenters = edges(1:end-1) + diff(edges)/2;
        plot(binCenters, counts, 'LineWidth', 1.5);
        hold on;
        xline(optParams(i), '--r', 'LineWidth', 1.5);
        title(paramNames(i));
        hold off;
    end
    if saveFigs
        saveas(gcf, "Figures/MCMC_hist_parallel_" + num2str(ntrials, "%.1e") + "trials.fig");
    end
    
    % Evolution of parameters
    figure(2);
    clf;
    hold on;
    sgtitle("Parallel MCMC: Parameter Evolution");
    for i = 1:length(paramNames)
        subplot(320 + i);
        plot(squeeze(models(i, randi(100,10,1), :)) ./ rescale_factor(i));
        yline(lb(i), 'LineWidth', 3);
        yline(ub(i), 'LineWidth', 3);
        ylim([lb(i), ub(i)]);
        title(paramNames(i));
    end
    if saveFigs
        saveas(gcf, "Figures/param_evolution_parallel_" + num2str(ntrials, "%.1e") + "trials.fig");
    end
    
    % Log likelihood trace
    figure(3);
    clf;
    plot(reshape(logP, 1, []));
    xlabel("Iterations")
    ylabel("Log Likelihood")
    if saveFigs
        saveas(gcf, "Figures/likelihood_parallel_" + num2str(ntrials, "%.1e") + "trials.fig");
    end
    
    % Top 10% correlation plot (example)
    num_top_elements = ceil(0.1 * length(reshape(logP, 1, [])));
    [~, sorted_indices] = sort(reshape(logP, 1, []), 'descend');
    top_indices = sorted_indices(1:num_top_elements);
    figure(4);
    clf;
    % Assuming: parameter 1 = HMM volume, 2 = dpHMM, 3 = vertical, 4 = horizontal, 6 = dpSC.
    vol_HMM = samples(:,1);
    dpHMM  = samples(:,2);
    vert   = samples(:,3);
    horiz  = samples(:,4);
    dpSC   = samples(:,6);
    plot(dpHMM .* vol_HMM, (4/3) * pi .* vert .* (horiz.^2) .* dpSC, '.')
    xlabel("dpHMM * HMM volume");
    ylabel("dpSC * SC volume");
    if saveFigs
        saveas(gcf, "Figures/volume_correl_parallel_" + num2str(ntrials, "%.1e") + "trials.fig");
    end

end

%% Local Log-Probability Function
function [lp, dummy] = local_log_prob(m, data, priormeans, N_iter, m_known, xopt, xtilt, yopt, ytilt, ...
                             zopt, invStdPWRL, tiltstd, tiltreduced, nanstatbeg, sigma, xbnds, rescale_factor)
    if any(m < xbnds(:,1)) || any(m > xbnds(:,2))
        lp = -Inf;
        dummy = 0;
        return;
    end
    % Call your model evaluation function.
    % Note: The order and number of additional arguments should match your original call.
    dprop = create_MCMC_data(m./rescale_factor, m_known, xopt, yopt, zopt, data, sigma, [xtilt, ytilt], nanstatbeg, ...
        priormeans, [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3), 1/tiltstd, 1/tiltstd]);
    % create_MCMC_data(m, m_guess, x, y, z, u, recstds, tilt, nanstatsbeg, priormeans, priorvariances)
    data_vec = data(:);    % Convert data to a column vector
    dprop_vec = dprop(:);  % Convert dprop to a column vector
    lp = -0.5 * norm((data_vec - dprop_vec) ./ sigma(:), 1) * 1e-2;
    dummy = 0;
end
