gps_weights1 = linspace(1,1e3, 10);
prior_weights1 = linspace(100, 2e3, 10);

gps_weights2 = linspace(1.5e2, 5e2, 10);
prior_weights2 = linspace(600, 1.5e3, 10);
% save("Data/l_curve_data_" + l_curve_type + "_dvinversion_5e5_2.mat", "l_curve_points", "l_curve_type", "dynamic_prior_weights", "dynamic_gps_weights", "posteriors_list");
% GPS L curve data
load("Data/l_curve_data_gps_dvinversion_5e5.mat", "l_curve_points", "l_curve_type", "posteriors_list");
l_curve_combined_gps = l_curve_points;
posteriors_combined_gps = posteriors_list;
load("Data/l_curve_data_gps_dvinversion_5e5_2.mat", "l_curve_points", "l_curve_type", "posteriors_list");
l_curve_combined_gps = [l_curve_combined_gps, l_curve_points];
posteriors_combined_gps = cat(3, posteriors_combined_gps, posteriors_list);

% Prior L curve data:
load("Data/l_curve_data_prior_dvinversion_5e5.mat", "l_curve_points", "l_curve_type", "posteriors_list");
l_curve_combined_prior = l_curve_points;
posteriors_combined_prior = posteriors_list;
load("Data/l_curve_data_prior_dvinversion_5e5_2.mat", "l_curve_points", "l_curve_type", "posteriors_list");
l_curve_combined_prior = [l_curve_combined_gps, l_curve_points];
posteriors_combined_prior = cat(3, posteriors_combined_gps, posteriors_list);

gps_weights_combined = [gps_weights1, gps_weights2];
prior_weights_combined = [prior_weights1, prior_weights2];


%%
% Set your ranges of interest here
target_gps_range = [100, 300];
target_prior_range = [800, 1200]; % Example: Prior weights between 800 and 1200

paramNames = ["dV HMM insar", "dV HMM GPS", "dV SC insar", "dV SC GPS"];

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
figure('Name', 'MCMC Posterior Histograms', 'Color', 'w', 'Position', [100, 100, 1200, 600]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Loop through the requested parameters (Param 1 and Param 2)
for i = 1:length(param_indices)
    p_idx = param_indices(i);
    nexttile;
    hold on;
    grid on;
    title("Parameter " + paramNames(i), 'FontSize', 14);
    xlabel('Parameter Value');
    ylabel('Frequency / Density');
    
    legend_entries = {};
    
    % Plot GPS-varying runs matching the criteria
    % We use a loop to plot each run individually so we can see the spread
    for k = 1:length(selected_gps_idx)
        run_idx = selected_gps_idx(k);
        data = posteriors_combined_gps(p_idx, :, run_idx);

        % Using histogram with transparency
        h = histogram(data, 50, 'Normalization', 'pdf', ...
            'DisplayStyle', 'bar', ...
            'FaceColor', [0 0.4470 0.7410], ... % Standard MATLAB Blue
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        if k == 1
            legend_entries{end+1} = 'GPS Vars'; % Add to legend only once
        else
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
    
    % Plot Prior-varying runs matching the criteria
    for k = 1:length(selected_prior_idx)
        run_idx = selected_prior_idx(k);
        data = posteriors_combined_prior(p_idx, :, run_idx);

        h = histogram(data, 50, 'Normalization', 'pdf', ...
            'DisplayStyle', 'bar', ...
            'FaceColor', [0.8500 0.3250 0.0980], ... % Standard MATLAB Red/Orange
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');

        if k == 1
            legend_entries{end+1} = 'Prior Vars';
        else
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
    
    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'best');
    end
    hold off;
end

sgtitle('Comparison of Selected MCMC Runs');