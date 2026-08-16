function plotHists(posterior, dp_posterior, subsamp_inds, opt_ind, optParams, lb, ub, paramNames, saveFigs)
%% Plot histogram of each parameter
% Toggle for delta V superscript labels: when true, label the InSAR-window
% and GPS-window volume changes as date range (1) and (2) instead of
% "InSAR"/"GPS". Applies to both HMM and SC delta V labels.
useDateRangeLabels = false;
% Smoothing knob for the ksdensity curves: >1 = smoother, <1 = more detail.
% Scales each parameter's automatically-chosen bandwidth (adapts to units).
smoothFactor = 1.8;
if useDateRangeLabels
    dvHMM_insar_sup = 'InSAR'; dvHMM_gps_sup = 'GPS';
    dvSC_insar_sup  = 'InSAR'; dvSC_gps_sup  = 'GPS';
else
    dvHMM_insar_sup = '\mathrm{InSAR}'; dvHMM_gps_sup = '\mathrm{GPS}';
    dvSC_insar_sup  = '\mathrm{InSAR}'; dvSC_gps_sup  = '\mathrm{GPS}';
end
plotParamNames = {
  ['$\Delta V_{\mathrm{HMM}}^{' dvHMM_insar_sup '}\ (\mathrm{km^3})$'], ...
  ['$\Delta V_{\mathrm{HMM}}^{' dvHMM_gps_sup '}\ (\mathrm{km^3})$'], ...
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
  ['$\Delta V_{\mathrm{SC}}^{' dvSC_insar_sup '}\ (\mathrm{km^3})$'], ...
  ['$\Delta V_{\mathrm{SC}}^{' dvSC_gps_sup '}\ (\mathrm{km^3})$'], ...
  '$\Delta P_{\mathrm{SC}}^{\mathrm{InSAR}}\ (\mathrm{MPa})$', ...
  '$\Delta P_{\mathrm{SC}}^{\mathrm{GPS}}\ (\mathrm{MPa})$', ...
  '$V_{\mathrm{SC}}\ (\mathrm{km}^3)$',
};
load("/Users/eric/Desktop/Stanford/Summer 23 Lab/Matlab/ForEric/no_vol_inversion/Data/paramDists.mat", "paramDists");
posterior_2 = posterior(:, subsamp_inds);
dp_posterior_2 = dp_posterior(:, subsamp_inds);
% Convert HMM top depth to centroid, include pressure changes:
taiyi_parameters = [1600.79, 914.47, 90, 0, 50, 200, -2.18e3-300, -4e7, ... 
     277.01, 1621.47, 63, 136, 0, 0, 0, -3630, -10e6];
ci = prctile(posterior(6,:), [5 95]) * 1e-3;
fprintf('%s: MAP = %.2f, 90%% CI [%.2f, %.2f]\n', 'top depth', 0.8, ci(1), ci(2));
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
optimizedM = get_full_m(zeros(1, 16), optParams', true, "insar");
optParams = [optParams(1:2); dp_posterior_2(1, opt_ind); dp_posterior_2(2, opt_ind); ...
    optParams(3:end-1); dp_posterior_2(end-2, opt_ind); dp_posterior_2(end-1, opt_ind); optParams(end)];
lb = [lb(1:2), -60e6, -60e6, lb(3:end-1), -60e6, -60e6, lb(end)];
ub = [ub(1:2), 0, 0, ub(3:end-1), 0, 0, ub(end)];
% Scale units to appropriate factor
unitScaling = [1e-9, 1e-9, 1e-6, 1e-6, 1e-9, 1e-3, 1e-3, 1e-3, 1, ...
    1e-3, 1e-3, 1e-3, 1, 1, 1, 1e-9, 1e-9, 1e-6, 1e-6, 1e-9];
histlims = [-8e-2, 0; -4e-2, 0; -40, 0; -40, 0; 1, 15; -0.1, 0.2; 0, 0.4; -3.5, -1.5; 1.3, 1.8; ...
    0, 2; -2.8, -2.0; -5, -3; 0.1, 0.6; 50, 80; 100, 180; -8e-2, 0; -4e-2, 0; -60, 0; -60, 0; 2, 15];
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

    if(i==8 || i == 12)
        mle = abs(mle);
        data = abs(data);
        tmpub = histlims(i, 1);
        histlims(i,1) = abs(histlims(i, 2));
        histlims(i,2) = abs(tmpub);
        posterior(i, :) = abs(posterior(i, :));
    end
    
    range_min = histlims(i,1); 
    range_max = histlims(i,2);
    xGrid = linspace(range_min, range_max, 200)';
    % Create Histogram
    s = unitScaling(i);
    [f_bg, xi_bg, bw_auto] = ksdensity(data(~isnan(data)) * s, xGrid);
    [f_bg, xi_bg] = ksdensity(data(~isnan(data)) * s, xGrid, 'Bandwidth', bw_auto * smoothFactor);
    post_peak = max(f_bg); % track peak posterior density for scaling the uniform box

    % Determine if current parameter is a post-processed pressure (indices 3, 4, 18, 19)
    is_derived = ismember(i, [3, 4, 18, 19]);
    if is_derived
        color_sub = [0.8500, 0.3250, 0.0980]; % MATLAB orange for derived
    else
        color_sub = [0, 0.4470, 0.7410];      % MATLAB blue for inverted
    end
        
    % Plot as a filled area
    area(xi_bg, f_bg, 'FaceColor', color_sub, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.15);
    hold on;
    % Plot Posterior
    hPost_sub = plot(xi_bg, f_bg, 'Color', color_sub, 'LineWidth', 3);
    
    % Save distinct handles for the legend
    if is_derived
        hPost_sub_der = hPost_sub;
    else
        hPost_sub_inv = hPost_sub;
    end
    
    % Plot second posterior
    if(~isempty(posterior))
        data2 = posterior(i, :);
        data2 = data2(~isnan(data2)) * s;
        valid_mask = (data2 >= range_min) & (data2 <= range_max);
        data2 = data2(valid_mask);
        [f_bg, xi_bg, bw_auto] = ksdensity(data2, xGrid);
        [f_bg, xi_bg] = ksdensity(data2, xGrid, 'Bandwidth', bw_auto * smoothFactor);
        post_peak = max(post_peak, max(f_bg));

        % Plot as a filled area
        hPost = area(xi_bg, f_bg, ...
            'FaceColor', [0.6 0.6 0.6], ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.3);
    end
    
    % Plot Priors
    % Extend prior beyond posterior range by 10% each side
    if(i == 8)
        lb_chosen = -3.5e3;
        ub_chosen = -1.5e3;
        priorName = "dcHMM";
        mle = abs(optimizedM(7)); % FIXED: added abs() here so it doesn't revert to negative
    end
    
    hMLE = xline(mle*s, '--r', 'LineWidth', 2);
    name = plotParamNames{i};
    
    % 1. Create the grid using the original true bounds (negative for depths)
    xGrid = linspace(lb_chosen, ub_chosen, 200)';

    dist_obj = paramDists.(priorName).dist;
    if strcmpi(dist_obj.DistributionName, 'Uniform')
        % Uniform prior: draw a dashed box spanning min to max of its bounds
        edges = [dist_obj.Lower, dist_obj.Upper];
        if(i == 8 || i == 12)
            edges = abs(edges);
        end
        lo = min(edges); hi = max(edges);
        % Scale the box height to sit just below the posterior peak so it is
        % visible on the same axis (rather than the true 1/(hi-lo) density).
        h = 0.9 * post_peak;
        % Clamp the box edges to the visible x-range so the vertical drops at
        % the limits stay on-screen (they get clipped otherwise, leaving only
        % the horizontal top edge visible). This renders a proper boxcar.
        lo_plot = max(lo*s, range_min);
        hi_plot = min(hi*s, range_max);
        hPrior = plot([lo_plot, lo_plot, hi_plot, hi_plot], [0, h, h, 0], '--', ...
            'Color', [0.3 0.3 0.3], 'LineWidth', 2.5);
    else
        % 2. Evaluate the PDF on the original negative grid
        p_prior = pdf(dist_obj, xGrid);

        % 3. Flip the x-coordinates to positive for plotting depths (indices 8 and 12)
        if(i == 8 || i == 12)
            xGrid = abs(xGrid);
        end

        % Rescale prior
        xGrid = xGrid * s;
        p_prior = p_prior / s;

        % Prior PDF in dark gray dashed
        if(~ismember(i, [1:5, 16:20]))
            hPrior = plot(xGrid, p_prior, '--', 'Color',[0.3 0.3 0.3], 'LineWidth',2.5);
        end
    end
   

    % Formatting
    grid off;
    set(gca, 'FontSize', 20, 'yticklabel', [], 'ytick', [], 'TickLength', [0.04, 0.04], 'LineWidth', 2);
    
    title(name, 'Interpreter','latex', 'FontSize', 30);
    axis square;
    xlim(histlims(i,:));
    
    % Print CI stats
    ci = prctile(data, [5 95]) * s;
    fprintf('%s: MAP = %.2f, 90%% CI [%.2f, %.2f]\n', name, mle*s, ci(1), ci(2));
end

% Add Legend with new handles for color distinction
lgd = legend([hPost, hPost_sub_inv, hPost_sub_der, hPrior, hMLE], ...
    {'Posterior PDF', 'Subsampled PDF (sampled)', 'Subsampled PDF (derived)', 'Prior PDF', 'MAP estimate'});
lgd.Layout.Tile = 'south'; % Moves legend to the bottom center shared space
lgd.Orientation = 'horizontal'; % Arranges items side-by-side
lgd.NumColumns = 3;
lgd.FontSize = 24;
lgd.Box = 'off';

% Add footnote to the bottom of the tiled layout explaining the asterisk
% xlabel(tl, '$^\ast$Derived from post-processing data', 'Interpreter', 'latex', 'FontSize', 22);

if saveFigs; exportgraphics(tl, './PaperFigs/combined_hist.jpg', 'Resolution', 300); end
end