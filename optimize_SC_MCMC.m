function [optParams, posterior, L_keep, gps_l2, insar_l2, prior_l2] = optimize_SC_MCMC(m_known, lb, ub, xopt, yopt, zopt, u1d, ...
     insarx, insary, insaru, look, insar_lengths, cinv_full, invStdPWRL, nanstatbeginning, ...
     ntrials, gps_weight, prior_weight, paramNames, burn, second_run, saveFigs)

GPS_std = 1./invStdPWRL;
priormeans = get_full_m(m_known, [], false, "insar");
if(second_run)
    priormeans = get_full_m(m_known, [], false, "insar", second_run);
    priormeans(1:2) = [4e9, 3e9]; 
end
out_of_bnds = priormeans > ub | priormeans < lb;
if(any(out_of_bnds)); error("Initial guess is not within bounds index " + strjoin(string(find(out_of_bnds)), ', ')); end
% priormeans = priormeans + 0.1*randn(1,6) .* priormeans;
% paramNames = ["HMM volume", "dpHMM", "vert semi-diameter", "horiz semi-diameter", "dip", "dpSC"];

bnds = [lb; ub]';
sigma = ones(size([u1d(:)])); % Insar stddev ~0.2m
% sigma = ones(size(insaru));
for i = 1:length(u1d(:))
    sigma(i) = GPS_std(ceil(i / size(u1d, 2)));
end

% [u1d(:);insaru]
% xstep = 0.02*ones(1,6); % 0.02
% xstep(4) = 0.007; % horiz semi-diam 

xstep = 1.4e-1*ones(1,size(bnds,1)); % 2.5e-2

[x_keep, L_keep, count, gps_l2, insar_l2, prior_l2] = mcmc('create_MCMC_data',[u1d(:);insaru(:)],priormeans,xstep, ...
    bnds, sigma, cinv_full, ntrials, gps_weight, prior_weight, paramNames, burn,...
    m_known, xopt, yopt, zopt, insarx, insary, look, insar_lengths, nanstatbeginning, second_run);

% burn = 4e3;

% Store results (SELECT LARGEST LL): 
[~, idx] = max(L_keep(burn:end));  % Index of the highest log-likelihood
posterior = x_keep(:, burn:end);
optParams = posterior(:, idx);
% if(subsample)
%     subsample_inds = [5, 6, 14, 15];
%     subsample_bnds = [0, 0.4e3; -2e3, -0.5e3; -10e6, 0; -10e6, 0];
%     subsample_lower = subsample_bnds(:,1);
%     subsample_upper = subsample_bnds(:,2);
%     posterior_subset = posterior(subsample_inds, :);
%     within_bounds = (posterior_subset >= subsample_lower) & (posterior_subset <= subsample_upper);
%     selection_mask = all(within_bounds, 1);
%     posterior = posterior(:,selection_mask);
% end

% Top 1% of loc likelihoods
% p = 0.01;                              % 5 % neighbourhood
% thresh = prctile(real(L_keep(burn:end)), 100*(1-p));
% idxTop = L_keep(burn:end) >= thresh;
% optParams = mean(real(posterior(:, idxTop)), 2);


%% Plotting

% figure(2);
% clf;
% hold on
% sgtitle("Step size = " + num2str(xstep));
% for i = 1:length(paramNames)
%     hold on;
%     subplot(4, 4, i); 
%     plot(x_keep(i,burn:end));
%     yline(lb(i), 'LineWidth', 3);
%     yline(ub(i), 'LineWidth', 3);
%     ylim([lb(i), ub(i)]);
%     title(paramNames(i));
%     hold off;
% end
% if(saveFigs); saveas(2, "./Figures/param_evolution_" + num2str(ntrials, "%.1e") + "trials.fig"); end
% 
% figure(3);
% clf;
% plot(real(L_keep));
% xlabel("Iterations")
% ylabel("Log Likelihood")
% if(saveFigs); saveas(3, "./Figures/likelihood_" + num2str(ntrials, "%.1e") + "trials.fig"); end


%% Get top 10% of estimates and plot the volume of HMM vs. volume of SC
% Find the number of elements corresponding to the top 10%
num_top_elements = ceil(0.1 * length(L_keep));

% Sort the list in descending order and get the indices
[~, sorted_indices] = sort(L_keep, 'descend');

% Get the indices of the top 10% values
top_indices = sorted_indices(1:num_top_elements);
L_keep = L_keep(burn:end);
% figure(4);
% clf;
% plot(x_keep(1, top_indices) .* x_keep(2, top_indices), (4/3) * pi * x_keep(3, top_indices) .* (x_keep(4, top_indices)).^2 .* x_keep(6, top_indices), '.')
% xlabel("dpHMM * vHMM");
% ylabel("dpSC * vSC");
% if(saveFigs); saveas(4, "Figures/volume_correl_" + num2str(ntrials, "%.1e") + "trials.fig"); end

end