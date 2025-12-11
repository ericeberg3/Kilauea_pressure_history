function plot_insar(insarx, insary, insaru, block_size, look, x, y, u1d, u1d_LSQ, nanstat, ...
    radscale, grid_res, GPSNameList, coast_new, cLimits, cmap, saveFigs)  
    %% InSAR Plot
        hold on;
        colormap(cmap);
        % Plot insar data
        currentColormap = colormap;
    
        % cLimits = [min(insaru), max(insaru)];
        nColors = size(currentColormap, 1);
        clim(cLimits);
        colorbar;
    
        % --- Overlay the quadtree squares ---
        % Determine grid resolution (assumes uniform spacing in insarx)
        grid_res = 31;  % 31 m/px resolution
        
        for i = 1:length(insarx)
            bs = block_size(i);       % block size (in grid cells)
            side_length = bs * grid_res; % convert block size to coordinate units
        
            % Compute lower-left corner from center (x(i), y(i))
            x_lower = insarx(i) - side_length/2;
            y_lower = insary(i) - side_length/2;
            
            % % Draw the square
            % rectangle('Position', [x_lower, y_lower, side_length, side_length], ...
            %           'EdgeColor', 'k', 'LineWidth', 1);
    
             % Map the displacement value to a color index
            colorIndex = round((insaru(i) - cLimits(1)) / (cLimits(2) - cLimits(1)) * (nColors - 1)) + 1;
            colorIndex = min(max(colorIndex, 1), nColors);  % ensure valid index
            faceColor = currentColormap(colorIndex, :);
        
            % Draw a filled square using patch. The vertices are given in clockwise order.
            patch([x_lower, x_lower+side_length, x_lower+side_length, x_lower], ...
              [y_lower, y_lower, y_lower+side_length, y_lower+side_length], ...
              faceColor, 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'LineWidth', 1, "HandleVisibility", "off");
        end
    
    
        % Get GPS equivalent insar displacements
        gpsInsarValues = u1d * look;
        colorIndices = round((gpsInsarValues - cLimits(1)) / (cLimits(2) - cLimits(1)) * (nColors - 1)) + 1;
        % Ensure indices are within bounds
        colorIndices = min(max(colorIndices, 1), nColors);
        gpsColors = currentColormap(colorIndices, :);
    
        % Plot GPS on top
        realquiver = quiver(x(~nanstat)', y(~nanstat)', u1d(:, 1) * radscale, u1d(:, 2) * radscale,  ...
            'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#4DBEEE', 'DisplayName', 'Data');
        % quiver3(mSCguess(5),mSCguess(6), mSCguess(7), optimizedM(13) - mSCguess(5), optimizedM(14) - mSCguess(6), optimizedM(15) - mSCguess(7), 'AutoScale', 'off', 'LineWidth', 2.75, 'MaxHeadSize', 0.5, 'Color', '#f77036', 'DisplayName','SC Center Shift');
        simquiver = quiver(x', y', u1d_LSQ(:, 1) * radscale, u1d_LSQ(:, 2) * radscale, ...
            'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#A2142F', 'DisplayName','Optimization Result');
        % scatter(x, y, 100, gpsColors, 'filled', 'HandleVisibility', 'off');
        scatter(x, y, 100, 'k', 'filled', 'HandleVisibility', 'off'); % Turned off color for now
    
        % Plot vert. displacement circles:
        theta = linspace(0, 2*pi, 50);  
        verticalScale = 0.3e4; 
        u1d_ind = 1;
        for k = 1:length(u1d_LSQ)
            if(~nanstat(k))
                line_style = '--';
                if(u1d(u1d_ind,3) < 0); line_style = '-'; end
    
                r = abs(u1d(u1d_ind,3)) * verticalScale;
                x_data = x(k) + r*cos(theta);
                y_data = y(k) + r*sin(theta);
    
                plot(x_data, y_data, line_style, 'Color', '#4DBEEE', 'LineWidth', 3, 'HandleVisibility','off');
                u1d_ind = u1d_ind + 1;
            end
    
            line_style = '--';
            if(u1d_LSQ(k,3) < 0); line_style = '-'; end
    
            r_lsq = abs(u1d_LSQ(k,3)) * verticalScale;
    
            x_sim = x(k) + r_lsq*cos(theta);
            y_sim = y(k) + r_lsq*sin(theta);
    
            plot(x_sim, y_sim, line_style, 'Color', '#A2142F', 'LineWidth', 2.5, 'HandleVisibility','off');
        end
    
        % Add reference circle for scale
        % Define reference location and scale
        x_ref = 6000;
        y_ref = 6000;
        r_ref = 0.1 * verticalScale;   % Reference vertical displacement magnitude (e.g., 1 unit)
        x_ref_circle = x_ref + r_ref*cos(theta);
        y_ref_circle = y_ref + r_ref*sin(theta);
        
        % Plot the reference circle
        plot(x_ref_circle, y_ref_circle, 'r', 'LineWidth', 2, 'HandleVisibility','off');
        text(x_ref - 270, y_ref - r_ref - 200, "10 cm", 'Color', 'red', 'FontSize', 12, 'HandleVisibility','off');
        quiver(x_ref, y_ref, -radscale * 0.5, 0, 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', 'red', 'HandleVisibility','off'); % Vector scale reference
        quiver(x_ref - 1e3, y_ref + 1e3, radscale * look(1), radscale*look(2), 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', 'red', 'HandleVisibility','off'); % Look vector direction
        
        text(x_ref - radscale*0.4, y_ref -radscale*0.1, "50 cm", 'Color', 'red', 'FontSize', 12, 'HandleVisibility','off');
        xlabel('Easting (m)', "FontSize", 18);
        ylabel('Northing (m)', "FontSize", 18);
        zlabel('Vertical (m)',  "FontSize", 18);
        title("Displacements at various stations (y oriented N, x oriented E)", "FontSize", 18)
        xlim([-8000, 8000]);
        ylim([-8000, 8000]);
        zlim([-8000, 8000]);
        GPSNameList(end + 1) = "SDH";
        
        cxy = llh2local(coast_new', [-155.2784, 19.4073]);
        cxy = cxy * 1000;
        plot(cxy(1, :)', cxy(2, :)', 'k.', 'HandleVisibility','off');
        
        % text((x-300)', (y+300)', (zeros(size(x)) + 400)', GPSNameList, 'FontSize', 16, 'HandleVisibility','off');
        % plot3(x', y', zeros(size(x))', '.', "Color", "#f76a05", "MarkerSize", 20, 'HandleVisibility','off')
        hold off;
        legend("FontSize", 18, "Location", "southwest");

end