function plotLcurve(l_curve_points, l_curve_type, prior_weights, gps_weights)
    %L curve plotting 
    % Prepare L-curve data, GPS vs. insar
    % xl = l_curve_points(1,:);
    % yl = l_curve_points(2,:);
    % xl = zeros(100, 1);
    % yl = zeros(100, 1);
    
    % L-curve data, GPS + insar vs. prior
    if(l_curve_type == "prior")
        xl = l_curve_points(1,:) + l_curve_points(2,:);
        yl = abs(l_curve_points(3,:));
        chosen_weights = prior_weights;
    else
        xl = l_curve_points(1,:);
        yl = l_curve_points(2,:);
        chosen_weights = gps_weights;
    end
    
    % Scatter, mapping color to gps_weights
    l_curve = figure(1); clf;
    
    % Choose colormap and add colorbar
    
    if(l_curve_type == "prior")
        L_scaling = (10 * 3 * 2e2 + 2379); % L_scaling = N_gps * w_gps + N_insar
        scatter_handle = scatter(xl(:), yl(:), 400, prior_weights(:), 'filled');
        % [~, ind] = min(abs(0.1 - prior_weights));
        colormap(parula);
        c = colorbar;
        c.Label.String = 'Prior Weight';
        xlabel("GPS L2 Norm + InSAR L2 Norm", 'FontSize', 24);
        ylabel("Prior Log Probability", 'FontSize', 24);
        % title("L-curve colored by prior weight", 'FontSize', 30);
    else
        scatter_handle = scatter(xl(:), yl(:), 400, gps_weights(:), 'filled');
        colormap(parula);
        c = colorbar;
        c.Label.String = 'GPS Weight';
        xlabel("GPS L2");
        ylabel("InSAR L2");
        % title("L-curve colored by GPS weight", 'FontSize', 30);
    end
    c.FontSize = 32;
    axis square;
    ax = gca;
    ax.FontSize = 32;
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [2 2 15 15]); 
    % xlim([420, 550])
    
    % Set up interactive data cursor
    dcm = datacursormode(gcf);
    set(dcm, 'UpdateFcn', @(~, event_obj) ...
        sprintf('GPS L2: %.2e\nInSAR L2: %.2e\nGPS Weight: %.1e', ...
                event_obj.Position(1), ...
                event_obj.Position(2), ...
                chosen_weights(event_obj.DataIndex)));
    
    exportgraphics(l_curve, "./PaperFigs/l_curve_" + l_curve_type + ".png", 'Resolution', 500);

end