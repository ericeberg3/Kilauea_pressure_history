
% save("Data/l_curve_data_" + l_curve_type + "_dvinversion_5e5_2.mat", "l_curve_points", "l_curve_type", "dynamic_prior_weights", "dynamic_gps_weights", "posteriors_list");
% GPS L curve data
load("../Data/l_curve_data_gps_dvinversion_5e5_nointersect.mat", "l_curve_points", "l_curve_type", "posteriors_list", "gps_weights");
gps_weights1 = gps_weights;
l_curve_combined_gps = l_curve_points;
posteriors_combined_gps = posteriors_list;
load("../Data/l_curve_data_gps_dvinversion_5e5_nointersect_2.mat", "l_curve_points", "l_curve_type", "posteriors_list", "gps_weights");
gps_weights2 = gps_weights;
l_curve_combined_gps = [l_curve_combined_gps, l_curve_points];
posteriors_combined_gps = cat(3, posteriors_combined_gps, posteriors_list);

% Prior L curve data:
load("../Data/l_curve_data_prior_dvinversion_5e5.mat", "l_curve_points", "l_curve_type", "posteriors_list", "prior_weights");
prior_weights1 = prior_weights;
l_curve_combined_prior = l_curve_points;
posteriors_combined_prior = posteriors_list;
load("../Data/l_curve_data_prior_dvinversion_5e5_2.mat", "l_curve_points", "l_curve_type", "posteriors_list", "prior_weights");
prior_weights2 = prior_weights;
l_curve_combined_prior = [l_curve_combined_gps, l_curve_points];
posteriors_combined_prior = cat(3, posteriors_combined_prior, posteriors_list);

gps_weights_combined = [gps_weights1, gps_weights2];
prior_weights_combined = [prior_weights1, prior_weights2];


%%
% Set your ranges of interest here
target_gps_range = [200 * 0.1, 260 * 1.9];
target_prior_range = [0.12 * 0.1, 0.12 * 1.9]; % Example: Prior weights between 800 and 1200

paramNames = ["\DeltaV^{InSAR}_{HMM}", "\DeltaV^{GPS}_{HMM}", "\DeltaV^{InSAR}_{SC}", "\DeltaV^{GPS}_{SC}"];

% Parameter indices to plot
param_indices = [1, 2, 14, 15]; 

% filtering
% Find indices in the weight arrays that satisfy the conditions
selected_gps_idx = find(gps_weights_combined >= target_gps_range(1) & ...
                        gps_weights_combined <= target_gps_range(2));

selected_prior_idx = find(prior_weights_combined >= target_prior_range(1) & ...
                          prior_weights_combined <= target_prior_range(2));

fprintf('Found %d runs in GPS range [%.1f, %.1f]\n', length(selected_gps_idx), target_gps_range);
fprintf('Found %d runs in Prior range [%.1f, %.1f]\n', length(selected_prior_idx), target_prior_range);

% Visualization
figure(1);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Loop through the requested parameters (Param 1 and Param 2)
for i = 1:length(param_indices)
    p_idx = param_indices(i);
    nexttile;
    hold on;
    title(paramNames(i), 'FontSize', 14);
    
    legend_entries = {};
    
    % Plot GPS-varying runs matching the criteria
    % We use a loop to plot each run individually so we can see the spread
    for k = 1:length(selected_gps_idx)
        run_idx = selected_gps_idx(k);
        data = posteriors_combined_gps(p_idx, :, run_idx)./1e9;

        % Using histogram with transparency
        h = histogram(data, 50, 'Normalization', 'pdf', ...
            'DisplayStyle', 'bar', ...
            'FaceColor', [0 0.4470 0.7410], ... % Standard MATLAB Blue
            'FaceAlpha', 0.2, ...
            'EdgeColor', 'none');

        if i == 1 && k == 1
            legend_entries{end+1} = 'GPS Variations'; % Add to legend only once
        else
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
    
    % Plot Prior-varying runs matching the criteria
    for k = 1:length(selected_prior_idx)
        run_idx = selected_prior_idx(k);
        data = posteriors_combined_prior(p_idx, :, run_idx)./1e9;

        h = histogram(data, 50, 'Normalization', 'pdf', ...
            'DisplayStyle', 'bar', ...
            'FaceColor', [0.8500 0.3250 0.0980], ... % Standard MATLAB Red/Orange
            'FaceAlpha', 0.4, ...
            'EdgeColor', 'none');

        if i == 1 && k == 1
            legend_entries{end+1} = 'Prior Variations';
        else
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
    set(gca, 'FontSize', 24)
    
    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'northwest');
    end
    hold off;
    box on;
    xlim([-inf, -0.005]);
end

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [2 2 21 12]); 
% sgtitle('Comparison of Selected MCMC Runs');
exportgraphics(gcf, "../PaperFigs/weight_perturbation.png", 'Resolution', 500);