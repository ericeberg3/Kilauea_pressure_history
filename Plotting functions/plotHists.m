function plotHists(posterior, optParams, lb, ub, paramNames, saveFigs)

%% Plot histogram of each parameter
plotParamNames = {
  '$\Delta V_{\mathrm{HMM}}^{\mathrm{InSAR}}\ (\mathrm{MPa})$', ...
  '$\Delta V_{\mathrm{HMM}}^{\mathrm{GPS}}\ (\mathrm{MPa})$', ...
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
  '$\Delta V_{\mathrm{SC}}^{\mathrm{InSAR}}\ (\mathrm{MPa})$', ...
  '$\Delta V_{\mathrm{SC}}^{\mathrm{GPS}}\ (\mathrm{MPa})$', ...
  '$V_{\mathrm{SC}}\ (\mathrm{km}^3)$',
};

load("/Users/eric/Desktop/Summer 23 Lab/Matlab/ForEric/no_vol_inversion/Data/paramDists.mat", "paramDists");

% Scale units to appropriate factor
unitScaling = [1e-9, 1e-9, 1e-9, 1e-3, 1e-3, 1e-3, 1, ...
    1e-3, 1e-3, 1e-3, 1, 1, 1, 1e-9, 1e-9, 1e-9];

histlims = [-50, 0; -40, 0; 0, 20; -0.1, 0.2; 0, 0.4; -1.5, -0.75; 1.3, 1.8; ...
    0, 2; -2.8, -2.0; -4.5, -3; 0.1, 0.6; 50, 80; 100, 180; -40, 0; -40, 0; 2, 15];

figure(5); clf;
tl = tiledlayout(4,4,'Padding','compact', 'TileSpacing','compact');

for i = 1:16
    ax = nexttile;
    
    % Retrieve data based on source map
    data = posterior(i, :); % From Calculated Pressures
    priorName = paramNames{i};
    mle = optParams(i);
    lb_chosen = lb(i);
    ub_chosen = ub(i);

    % Create Histogram
    s = unitScaling(i);
    [counts, edges] = histcounts(data(~isnan(data)), 50, 'Normalization', 'pdf');
    binCenters = (edges(1:end-1) + diff(edges)/2) * s;
    counts = counts / s;
    
    % Plot Posterior
    hPost = plot(binCenters, counts, 'Color',[0 0.4470 0.7410], 'LineWidth', 3);
    hold on;
    
    % Plot MLE / Mean (For pressure/vol, use mean of samples; for geom use optParams)
    val_mle = optParams(i) * s;
    xline(mle*s, '--r', 'LineWidth', 2);
    
    % Plot Priors
    % Extend prior beyond posterior range by 10% each side
    span = edges(end) - edges(1);
    xGrid = linspace(lb_chosen, ub_chosen, 400)';
    name = plotParamNames{i};
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
    title(name, 'Interpreter','latex', 'FontSize', 16);
    axis square;
    
    % Print CI stats
    ci = prctile(data, [5 95]) * s;
    fprintf('%s: 90%% CI [%.2f, %.2f]\n', name, ci(1), ci(2));
end

% Add Legend
if saveFigs; exportgraphics(tl, './PaperFigs/combined_hist.png', 'Resolution', 300); end

end