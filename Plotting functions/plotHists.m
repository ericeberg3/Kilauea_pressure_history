function plotHists(posterior, dp_posterior, subsamp_inds, opt_ind, optParams, lb, ub, paramNames, saveFigs)

%% Plot histogram of each parameter
plotParamNames = {
  '$\Delta V_{\mathrm{HMM}}^{\mathrm{InSAR}}\ (\mathrm{km^3})$', ...
  '$\Delta V_{\mathrm{HMM}}^{\mathrm{GPS}}\ (\mathrm{km^3})$', ...
  '$\Delta P_{\mathrm{HMM}}^{\mathrm{InSAR}}\ (\mathrm{MPa})$', ...
  '$\Delta P_{\mathrm{HMM}}^{\mathrm{GPS}}\ (\mathrm{MPa})$', ...
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
  '$\Delta V_{\mathrm{SC}}^{\mathrm{InSAR}}\ (\mathrm{km^3})$', ...
  '$\Delta V_{\mathrm{SC}}^{\mathrm{GPS}}\ (\mathrm{km^3})$', ...
  '$\Delta P_{\mathrm{SC}}^{\mathrm{InSAR}}\ (\mathrm{MPa})$', ...
  '$\Delta P_{\mathrm{SC}}^{\mathrm{GPS}}\ (\mathrm{MPa})$', ...
  '$V_{\mathrm{SC}}\ (\mathrm{km}^3)$',
};

load("/Users/eric/Desktop/Summer 23 Lab/Matlab/ForEric/no_vol_inversion/Data/paramDists.mat", "paramDists");
posterior_2 = posterior(:, subsamp_inds);
dp_posterior_2 = dp_posterior(:, subsamp_inds);

% Convert HMM top depth to centroid, include pressure changes:
taiyi_parameters = [1600.79, 914.47, 90, 0, 50, 200, -2.18e3-300, -4e7, ... 
     277.01, 1621.47, 63, 136, 0, 0, 0, -3630, -10e6];
for i = 1:length(posterior)
    m_tmp = get_full_m(taiyi_parameters, posterior(:,i)', true, "insar");
    posterior(6,i) = m_tmp(7);
end
for i = 1:length(posterior_2)
    m_tmp = get_full_m(taiyi_parameters, posterior_2(:,i)', true, "insar");
    posterior_2(6,i) = m_tmp(7);
end

% Add dp terms to the posterior and posterior_2
posterior = [posterior(1:2, :); dp_posterior(1:2, :); posterior(3:end-1, :); ...
    dp_posterior(end-2:end-1, :); posterior(end,:)];
posterior_2 = [posterior_2(1:2, :); dp_posterior_2(1:2, :); posterior_2(3:end-1, :); ...
    dp_posterior_2(end-2:end-1, :); posterior_2(end,:)];
paramNames_full = [paramNames{1:2}, "dpHMM_insar", "dpHMM_gps", paramNames{3:end-1}, ...
    "dpSC_insar", "dpSC_gps", paramNames{end}];
optParams = [optParams(1:2); dp_posterior_2(1, opt_ind); dp_posterior_2(2, opt_ind); ...
    optParams(3:end-1); dp_posterior_2(end-2, opt_ind); dp_posterior_2(end-1, opt_ind); optParams(end)];
lb = [lb(1:2), -60e6, -60e6, lb(3:end-1), -60e6, -60e6, lb(end)];
ub = [ub(1:2), 0, 0, ub(3:end-1), 0, 0, ub(end)];

% Scale units to appropriate factor
unitScaling = [1e-9, 1e-9, 1e-6, 1e-6, 1e-9, 1e-3, 1e-3, 1e-3, 1, ...
    1e-3, 1e-3, 1e-3, 1, 1, 1, 1e-9, 1e-9, 1e-6, 1e-6, 1e-9];

histlims = [-8e-2, 0; -4e-2, 0; -40, 0; -40, 0; 0, 20; -0.1, 0.2; 0, 0.4; -3.5, -1.5; 1.3, 1.8; ...
    0, 2; -2.8, -2.0; -4.5, -3; 0.1, 0.6; 50, 80; 100, 180; -8e-2, 0; -4e-2, 0; -60, 0; -60, 0; 2, 15];

figure(5); clf;
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 2000, 1600]);
tl = tiledlayout(4,5,'Padding','compact', 'TileSpacing','compact');

for i = 1:20
    ax = nexttile;
    
    % Retrieve data based on source map
    data = posterior_2(i, :); % From Calculated Pressures
    priorName = paramNames_full(i);
    mle = optParams(i);
    lb_chosen = lb(i);
    ub_chosen = ub(i);
    
    range_min = histlims(i,1); 
    range_max = histlims(i,2);
    xGrid = linspace(range_min, range_max, 200)';
    % Create Histogram
    s = unitScaling(i);
    % [counts, edges] = histcounts(data(~isnan(data)), 30, 'Normalization', 'pdf');
    [f_bg, xi_bg] = ksdensity(data(~isnan(data)) * s, xGrid);
        
    % Plot as a filled area
    area(xi_bg, f_bg, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.15);
    hold on;
    % Plot Posterior
    hPost_sub = plot(xi_bg, f_bg, 'Color',[0 0.4470 0.7410], 'LineWidth', 3);
    
    % histogram(data(~isnan(data)) * s, 50, 'Normalization', 'pdf', ...
    %     'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    

    % Plot second posterior
    if(~isempty(posterior))
        data2 = posterior(i, :);
        data2 = data2(~isnan(data2)) * s;

        % hPost = histogram(data2, 50, 'Normalization', 'pdf', ...
        %     'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        valid_mask = (data2 >= range_min) & (data2 <= range_max);
        data2 = data2(valid_mask);
        [f_bg, xi_bg] = ksdensity(data2, xGrid);
        
        % Plot as a filled area
        hPost = area(xi_bg, f_bg, ...
            'FaceColor', [0.6 0.6 0.6], ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.3);
    end
    
    % Plot MLE / Mean (For pressure/vol, use mean of samples; for geom use optParams)
    hMLE = xline(mle*s, '--r', 'LineWidth', 2);
    
    % Plot Priors
    % Extend prior beyond posterior range by 10% each side
    if(i == 8)
        lb_chosen = -3.5e3;
        ub_chosen = -1.5e3;
        priorName = "dcHMM";
    end
    % span = edges(end) - edges(1);
    name = plotParamNames{i};
    xGrid = linspace(lb_chosen, ub_chosen, 200)';
    p_prior = pdf(paramDists.(priorName).dist, xGrid);
    
    % Rescale prior 
    xGrid = xGrid * s;
    p_prior = p_prior / s;
    % p_prior = max(data(~isnan(data)) * s) * p_prior/(max(p_prior));
    
    % Prior PDF in dark gray dashed 
    if(~ismember(i, [1:4, 16:20]))
        hPrior = plot(xGrid, p_prior, '--', 'Color',[0.3 0.3 0.3], 'LineWidth',2.5);
    end

    % xlim([lb_chosen * s, ub_chosen * s]);
    % Formatting
    grid off;
    set(gca, 'FontSize', 16);
    title(name, 'Interpreter','latex', 'FontSize', 24);
    axis square;
    xlim(histlims(i,:));
    
    % Print CI stats
    ci = prctile(data, [5 95]) * s;
    fprintf('%s: 90%% CI [%.2f, %.2f]\n', name, ci(1), ci(2));
end

% Add Legend
lgd = legend([hPost, hPost_sub, hPrior, hMLE], {'Posterior PDF', 'Subsampled PDF', 'Prior PDF', 'MAP estimate'});
lgd.Layout.Tile = 'south'; % Moves legend to the bottom center shared space
lgd.Orientation = 'horizontal'; % Arranges items side-by-side
lgd.FontSize = 24;
lgd.Box = 'off';

if saveFigs; exportgraphics(tl, './PaperFigs/combined_hist.png', 'Resolution', 300); end

end