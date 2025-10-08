function makeplots(x, y, GPS_llh, u, u1d, ux, uy, uz, u_low, u_high, insarx, insary, insaru, insaru_pred, block_size, look, tiltx, tilty, ...
    usim, t, nanstat, nanstatbeginning, finalindex, collapset, ...
    dp, dp_low, dp_high, tau, optParams, optimizedM, GPSNameList, gTiltHMM, gTiltSC, xtilt, ytilt, tiltreduced, radscale, ...
    coast_new, taiyi_parameters, disptype, ntrials, offsets, saveFigs)

    % u1d = squeeze(u(:, :, end-finalindex));
    % nanstat = isnan(u1d(:, 1));
    % u1d = u1d(~nanstat, :);
    %% Plotting with solved offsets
    ux = ux(:, :) + offsets(:, 1);
    uy = uy(:, :) + offsets(:, 2);
    uz = uz(:, :) + offsets(:, 3);

    %% Plots
    % Convert time into matlab dateyear
    year = floor(t);
    partialYear = mod(t,1);
    date0 = datenum(num2str(year),'yyyy');
    date1 = datenum(num2str(year+1),'yyyy');
    daysInYear = date1 - date0;
    t = date0 + partialYear .* daysInYear;
    t = datetime(t, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy');

    year = floor(collapset);
    partialYear = mod(collapset,1);
    date0 = datenum(num2str(year),'yyyy');
    date1 = datenum(num2str(year+1),'yyyy');
    daysInYear = date1 - date0;
    collapset = date0 + partialYear .* daysInYear;
    collapset = datetime(collapset, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy');
    
    % %% Pressure vs. time figure with dual y-axes
    % fig = figure(10); clf; 
    % 
    % % Define scale factor and offset for SC pressure (adjust these as needed)
    % sc_scale = 0.5;   % Scaling factor for SC data
    % sc_offset = 5;  % Downward offset for SC data
    % 
    % hmm_scale = 0.5;
    % 
    % % Plot HMM pressure on the right y-axis
    % ax1 = axes;  % This axis will host the HMM data and the transformed SC data
    % hold(ax1, 'on');
    % set(ax1, 'YColor', '#0072BD')
    % 
    % hmm_pressure = dp(1:end-finalindex,1) * optimizedM(8)/1e6;
    % hmm_pressure_adj = hmm_scale * hmm_pressure;
    % h1 = plot(ax1, t(1:end-finalindex), hmm_pressure_adj, ...
    %     'Color', '#0072BD', 'DisplayName', 'HMM Pressure', 'LineWidth', 4);
    % ylabel('HMM Pressure (MPa)');
    % hold on;
    % 
    % % Plot HMM confidence intervals on the right axis
    % xPatch = [t(1:end-finalindex); flipud(t(1:end-finalindex))];
    % yPatch_HMM = [hmm_scale * dp_low(1:end-finalindex, 1); flipud(hmm_scale * dp_high(1:end-finalindex, 1))];
    % 
    % patch(ax1, xPatch, yPatch_HMM, 'blue', 'FaceAlpha', 0.15, 'EdgeColor', 'none'); % HMM conf interval
    % 
    % % Adjust the SC pressure data with scaling and offset
    % sc_pressure = dp(1:end-finalindex, 2)*optimizedM(16)/(1e6);
    % sc_pressure_adj = sc_scale * sc_pressure - sc_offset;
    % h2 = plot(ax1, t(1:end-finalindex), sc_pressure_adj, 'DisplayName', 'SC Pressure', 'LineWidth', 4);
    % 
    % % Plot SC confidence intervals on the left axis
    % yPatch_SC = [sc_scale*(dp_low(1:end-finalindex, 2)) - sc_offset; ...
    %              flipud(sc_scale*(dp_high(1:end-finalindex, 2)) - sc_offset)];
    % patch(ax1, xPatch, yPatch_SC, 'red', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    % ylim([-20, 5]);
    % 
    % % Create axes label for SC
    % ax1_pos = get(ax1, 'Position');
    % ax1_pos(1) = ax1_pos(1) - 0.05;
    % ax2 = axes('Position', ax1_pos, ...
    %        'Color', 'none', ...
    %        'YAxisLocation', 'left', ...
    %        'XAxisLocation', 'top', ... 
    %        'XTick', [], ...       
    %        'Box', 'off');
    % 
    % set(ax2, 'YLim', [-20, 5], 'YColor','#D95319');
    % ax2.XAxis.Visible = 'off';
    % 
    % ticks = get(ax2, 'YTick');
    % tickLabels = (ticks + sc_offset) / sc_scale;
    % set(ax2, 'YTickLabel', tickLabels);
    % ylabel(ax2, 'SC Pressure (MPa)');
    % 
    % ticks = get(ax1, 'YTick');
    % tickLabels = ticks/hmm_scale;
    % set(ax1, 'YTickLabel', tickLabels);
    % 
    % % Set common figure properties
    % set(ax2, 'FontSize', 20);
    % set(ax1, 'FontSize', 20);
    % % title("Pressure vs. Time");
    % xlabel('Time');
    % legend([h1, h2], 'Location', 'northeast', "FontSize", 24);
    % hold off
    % 
    % 
    % disp("HMM Total Pressure Drop: " + dp(end-finalindex, 1)*optimizedM(8)/(1e6) + " ub: " + dp_high(end-finalindex, 1) + " lb: " + dp_low(end-finalindex,1));
    % disp("SC Total Pressure Drop: " + dp(end-finalindex, 2)*optimizedM(16)/(1e6) + " ub: " + dp_high(end-finalindex, 2) + " lb: " + dp_low(end-finalindex,2));
    % if(saveFigs)
    %     % saveas(10, "./Figures/pvt_" + num2str(ntrials, "%.1e") + "trials.fig");
    %     exportgraphics(fig, './PaperFigs/pvt.png', 'Resolution', 500);
    % end
    %% Subplot pvt version
    
    fig = figure(11);
    
    % --- global text size & figure size --------------------------------------
    baseFont = 22;                           % main tick‑label size
    % set(fig,'Color','w',...
    %         'Units','centimeters','Position',[4 4 20 10],...   % a bit wider
    %         'DefaultAxesFontSize', baseFont,...
    %         'DefaultAxesTitleFontSizeMultiplier', 1.25,...
    %         'DefaultAxesLabelFontSizeMultiplier', 1.25,...
    %         'DefaultLegendFontSize', baseFont);
    
    % --------- keep the original scaling variables but neutralise them -------
    hmm_scale = 1;
    sc_scale  = 1;
    sc_offset = 0;
    
    idx  = 1:numel(t)-finalindex;
    tvec = t(idx);
    
    % Pre‑compute pressures in MPa
    hmm_pressure = dp(idx,1) * optimizedM(8)  / 1e6;
    sc_pressure  = dp(idx,2) * optimizedM(16) / 1e6;
    
    % Layout
    tlo = tiledlayout(1,3,'TileSpacing','tight','Padding','tight');
    
    % ======================== 1) HMM subplot =================================
    ax1 = nexttile;  hold(ax1,'on'); 
    
    % <<<<<< ORIGINAL PATCH / ERROR‑BAR CODE >>>>>>
    xPatch     = [tvec; flipud(tvec)];
    yPatch_HMM = [hmm_scale * dp_low(idx,1); ...
                  flipud(hmm_scale * dp_high(idx,1))];
    greyColor = [0.7, 0.7, 0.7];
    pressure_color = "#c70000";
    patch(ax1, xPatch, yPatch_HMM, greyColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    % <<<<<< END COPY >>>>>>
    
    plot(ax1, tvec, hmm_scale*hmm_pressure,...
         'Color',pressure_color, 'LineWidth',2,'DisplayName','HMM Pressure');
    ax1.XAxis.FontSize = baseFont - 4;
    ax1.YAxis.FontSize = baseFont - 4;
    ylabel('HMM Pressure (MPa)','FontSize',baseFont+2)
    title('Halemaʻumaʻu (HMM)', 'FontSize', baseFont+4)
    grid on; box on; axis square
    set(ax1,'LineWidth',1.4, 'YLim', [-15, 5]);

    % ======================== 2) SC subplot ==================================
    ax2 = nexttile;  hold(ax2,'on'); 
    
    yPatch_SC = [sc_scale * dp_low(idx,2) - sc_offset; ...
                 flipud(sc_scale * dp_high(idx,2) - sc_offset)];
    patch(ax2, xPatch, yPatch_SC, greyColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
    plot(ax2, tvec, sc_scale*sc_pressure - sc_offset,...
         'Color',pressure_color,'LineWidth',2,'DisplayName','SC Pressure');
    ax2.XAxis.FontSize = baseFont - 4;
    ax2.YAxis.FontSize = baseFont - 4;
    ylabel('SC Pressure (MPa)','FontSize',baseFont+2)
    title('South Caldera (SC)', 'FontSize', baseFont+4)
    grid on; box on; axis square;
    set(ax2,'LineWidth',1.4)


    % ======================== 3) tau subplot ==================================
    ax3 = nexttile; hold(ax3, 'on');
    plot(ax3, tvec, tau(idx) .* 1e-6, 'LineWidth', 2);
    grid on; box on; axis square;
    set(ax3, 'LineWidth', 1.4);
    ax3.XAxis.FontSize = baseFont - 4;
    ax3.YAxis.FontSize = baseFont - 4;
    ylabel("Average Shear Stress (MPa)", 'FontSize', baseFont+2);
    title("Average shear stress vs. time", 'FontSize', baseFont+4)
    
    %%  Collapse markers & amplitude calculation
    %   Put this block *after* the two plots are finished (but before exportgraphics).
    
    delta_t   = hours(1.5);
    colStart  = [0.4 0.4 0.4];  % dark‑grey   – collapse start
    colSample = [0.7 0.7 0.7];  % light‑grey  – sample point
    
    ampl_HMM = nan(numel(collapset),1);
    ampl_SC = nan(numel(collapset),1);
    ampl_tau = nan(numel(collapset),1);
    
    Z_SCORE_90_10 = 1.2816;
    for k = 1:numel(collapset)
        t0 = collapset(k) - hours(1);
        t1 = t0 + delta_t;            % time at which we sample the “post‑collapse” value
        
        % -------- vertical lines on both subplots ---------------------------
        % collapse start
        % xline(ax1, t0, '--', 'Color', colStart,  'LineWidth', 1.3);
        % xline(ax2, t0, '--', 'Color', colStart,  'LineWidth', 1.3);
        % % sample point
        % xline(ax1, t1, ':',  'Color', colSample, 'LineWidth', 1.3);
        % xline(ax2, t1, ':',  'Color', colSample, 'LineWidth', 1.3);
        
        % -------- amplitude (HMM) ------------------------------------------
        % find nearest indices in time vector
        [~, i0] = min(abs(tvec - t0));
        [~, i1] = min(abs(tvec - t1));
        
        if i1 <= numel(tvec)          % make sure the sample point exists
            ampl_HMM(k, 1) = hmm_pressure(i1) - hmm_pressure(i0);
            ampl_SC(k, 1) = sc_pressure(i1) - sc_pressure(i0);
            ampl_tau(k,1) = abs(tau(i1) - tau(i0));
        end
    end
    
    avgAmpl_HMM = median(abs(ampl_HMM),1, 'omitnan');
    avgAmpl_SC = median(abs(ampl_SC), 1, 'omitnan');
    ampl_ratio = abs(ampl_SC(:,1) ./ ampl_HMM(:,1));
    fprintf('Average HMM collapse amplitude  %.3f MPa, pressure drop = %.3f MPa, [%.3f, %.3f] \n', ...
        avgAmpl_HMM(1), hmm_pressure(end), dp_low(end,1), dp_high(end,1));
    fprintf('Average SC collapse amplitude  %.3f MPa, pressure drop = %.3f MPa, [%.3f, %.3f] \n', ...
        avgAmpl_SC(1), sc_pressure(end), dp_low(end,2), dp_high(end,2));
    fprintf('Amplitude ratio (SC/HMM) = %.3f \n', mean(ampl_ratio, 'omitnan'))
    fprintf(' SC pressure drop = %.3f MPa, [%.3f, %.3f] \n', sc_pressure(end), ...
        dp_low(end,2), dp_high(end,2));
    fprintf('Average shear stress drop: %.2e \n', median(ampl_tau, 'omitnan'));
    
    
    % ------------------------ Export (optional) ------------------------------
    if saveFigs
        exportgraphics(fig,'./PaperFigs/pressure_histories_split.png','Resolution',600)
    end


    %% Making grid of displacements and tilt
    % for j = 1:1
    %     disptype = j;
    %     fig = figure(7);
    %     % disptype = 2; % 1 = x, 2 = y, 3 = z
    % 
    %     tlo = tiledlayout(4,4);
    %     if(disptype == 1)
    %         title(tlo, "East Displacement vs. Time", 'FontSize', 24);
    %     elseif(disptype == 2)
    %         title(tlo, "North Displacement vs. Time", 'FontSize', 24);
    %     else
    %         title(tlo, "Vertical Displacement vs. Time", 'FontSize', 24);
    %     end
    % 
    %     for i = 1:max(size(GPSNameList) + 1)
    %         nexttile
    %         if(i < length(GPSNameList) + 1)
    %             if(disptype == 1)
    %                 plot(t(1:end-finalindex), ux(i, 1:end-finalindex), '-', 'DisplayName', 'GPS', 'LineWidth', 1.2);
    %                 hold on;
    %                 plot(t(1:end-finalindex), usim(1:end-finalindex, 1, i), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
    %                 plot(t(end - finalindex), ux(i, end - finalindex), 'o-', 'MarkerFaceColor','red');
    % 
    %                 plot([0, 0], [0, offsets(i, j)], "LineWidth", 8);
    %             elseif(disptype == 2)
    %                 plot(t(1:end-finalindex), uy(i, 1:end-finalindex), '-', 'DisplayName', 'GPS', 'LineWidth', 1.2);
    %                 hold on;
    %                 plot(t(1:end-finalindex), usim(1:end-finalindex, 2, i), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
    %                 plot(t(end - finalindex), uy(i, end - finalindex), 'o-', 'MarkerFaceColor','red');
    % 
    %                 plot([0, 0], [0, offsets(i, j)], "LineWidth", 8);
    %             else
    %                 plot(t(1:end-finalindex), uz(i, 1:end-finalindex), '-', 'DisplayName', 'GPS', 'LineWidth', 1.2);
    %                 hold on;
    %                 plot(t(1:end-finalindex), squeeze(usim(1:end-finalindex, 3, i)), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
    %                 plot(t(end - finalindex), uz(i, end - finalindex), 'o-', 'MarkerFaceColor','red');
    % 
    %                 plot([0, 0], [0, offsets(i, j)], "LineWidth", 8);
    %                 % xline(collapset);
    %             end
    %             if(GPSNameList(i) == "UWEV" || GPSNameList(i) == "BYRL" || GPSNameList(i) == "CRIM")
    %                     set(gca,'Color','k');
    %             end
    %             ylabel("Displacement (m)");
    %             title(GPSNameList(i))
    %         else
    %             simtiltx = (gTiltHMM(1) .* dp(:, 1)) + (gTiltSC(1) .* dp(:, 2));
    %             plot(t(1:end-finalindex),tiltx(1:end-finalindex), '-', 'DisplayName', 'Data', 'LineWidth', 1.2);
    %             hold on;
    %             plot(t(1:end-finalindex), simtiltx(1:end-finalindex), '-', 'DisplayName', 'LSQ', 'LineWidth', 1.6);
    %             title("Tilt e");
    %             ylim([-30, 250]);
    %             ylabel("Tilt (µrad)")
    %             hold off;
    %             nexttile;
    %             simtilty = (gTiltHMM(2) .* dp(:, 1)) + (gTiltSC(2) .* dp(:, 2));
    %             plot(t(1:end-finalindex), tilty(1:end-finalindex), '-', 'DisplayName', 'Data', 'LineWidth', 1.2);
    %             hold on;
    %             plot(t(1:end-finalindex), simtilty(1:end-finalindex), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
    %             ylim([-90, 130]);
    %             ylabel("Tilt (µrad)")
    %             title("Tilt n");
    %         end
    %         hold off;
    %     end
    %     leg = legend('Orientation', 'Horizontal');
    %     leg.Layout.Tile = 'north';
    %     leg.FontSize = 14;
    %     if(saveFigs)
    %         % saveas(7, "./Figures/displacements_" + disptype + "_" + num2str(ntrials, "%.1e") + "trials.fig"); 
    %         exportgraphics(fig, "./PaperFigs/displacements_" + disptype + ".png", 'Resolution', 500);
    %     end
    % end


% —— Toggle inclusion of the three inactive stations ——
includeInactive = true;  % set to true to include UWEV, BYRL, CRIM

% —— Identify which stations to plot ——
inactive = ismember(GPSNameList,["UWEV","BYRL","CRIM"]);
allIdx   = 1:numel(GPSNameList);
if includeInactive
    plotIdx = allIdx;
else
    plotIdx = allIdx(~inactive);
end
nPlot    = numel(plotIdx);

% —— SDH (tilt) row? ——
hasSDH = true;                 % we’ll add it explicitly
nRows  = nPlot + hasSDH;       % add a row for SDH

% —— Shared settings ——
comps      = {'East','North','Up'};
ts         = t(1:end-finalindex);         % common time vector
dataCol    = [0 0.4470 0.7410];           % blue
fitCol     = [0.8500 0.3250 0.0980];      % red

% —— Precompute y‐limits per GPS station (unchanged) —— 
yLims = zeros(nPlot,2);
for row = 1:nPlot
    iSite = plotIdx(row);

    lo =  inf; hi = -inf;
    for col = 1:3
        low   = squeeze(u_low(iSite, col, 1:end-finalindex));
        high  = squeeze(u_high(iSite, col, 1:end-finalindex));
        switch col
            case 1
                dtemp = ux(iSite,1:end-finalindex);
                mtemp = squeeze(usim(1:end-finalindex,1,iSite));
            case 2
                dtemp = uy(iSite,1:end-finalindex);
                mtemp = squeeze(usim(1:end-finalindex,2,iSite));
            case 3
                dtemp = uz(iSite,1:end-finalindex);
                mtemp = squeeze(usim(1:end-finalindex,3,iSite));
        end
        lo = min([lo; low(:); dtemp(:); mtemp(:)]);
        hi = max([hi; high(:); dtemp(:); mtemp(:)]);
    end
    pad = (hi - lo) * 0.05;  % 5% padding
    yLims(row,:) = [lo - pad, hi + pad];
end

% —— SDH y-limits from tiltx/tilty and last-comp band —— 
if hasSDH
    iTilt = size(u_low,1);  % last station index holds SDH band
    lowSDH  = squeeze(u_low(end, 3, 1:end-finalindex));
    highSDH = squeeze(u_high(end, 3, 1:end-finalindex));

    dE = tiltx(1:end-finalindex);
    dN = tilty(1:end-finalindex);

    % If model predictions exist at the same last index, plot them; else NaN
    if ndims(usim) >= 3 && size(usim,3) >= iTilt
        mE = squeeze(usim(1:end-finalindex,1,iTilt));
        mN = squeeze(usim(1:end-finalindex,2,iTilt));
    else
        mE = nan(size(ts));
        mN = nan(size(ts));
    end

    lo = min([lowSDH(:); dE(:); mE(:); dN(:); mN(:)]);
    hi = max([highSDH(:); dE(:); mE(:); dN(:); mN(:)]);
    pad = (hi - lo) * 0.05;
    yLimsSDH = [lo - pad, hi + pad];
end

% —— Create compact grid ——
fig = figure(8); clf;
set(fig, 'Units','normalized','Position',[0.1 0.1 0.6 0.8]);
tl = tiledlayout(nRows,3, 'TileSpacing','none','Padding','none');

% —— Plot all regular GPS stations (unchanged behavior) —— 
for row = 1:nPlot
    iSite = plotIdx(row);

    % Check if the current station is inactive
    isInactive = ismember(GPSNameList(iSite), ["UWEV","BYRL","CRIM"]);
    
    % Choose the color based on station status
    if isInactive
        label_color = 'r';
    else
        label_color = 'k';
    end

    for col = 1:3
        ax = nexttile((row-1)*3 + col);
        hold(ax,'on');
        % Alternate y-axis side per tile
        if mod((row-1)*3 + col, 2) == 0
            ax.YAxisLocation = 'right';
        else
            ax.YAxisLocation = 'left';
        end

        low   = squeeze(u_low(iSite, col, 1:end-finalindex));
        high  = squeeze(u_high(iSite, col, 1:end-finalindex));
        ts_row  = ts(:)'; low_row = low(:)'; hi_row = high(:)';
        x_fill = [ts_row, fliplr(ts_row)];
        y_fill = [low_row, fliplr(hi_row)];
        fill(ax, x_fill, y_fill, fitCol, ...
             'EdgeColor','none', 'FaceAlpha',0.6, 'HandleVisibility','off');

        switch col
            case 1
                d = ux(iSite,1:end-finalindex);
                m = squeeze(usim(1:end-finalindex,1,iSite));
            case 2
                d = uy(iSite,1:end-finalindex);
                m = squeeze(usim(1:end-finalindex,2,iSite));
            case 3
                d = uz(iSite,1:end-finalindex);
                m = squeeze(usim(1:end-finalindex,3,iSite));
        end
        hdata = plot(ax, ts, d, '-', 'Color', dataCol, 'LineWidth',1.2);
        hfit  = plot(ax, ts, m, '-', 'Color', fitCol,  'LineWidth',1.6);

        xline(ax, 0, 'k:','LineWidth',1);
        grid(ax,'on'); ax.FontSize = 16;

        ax.XTickLabel = [];

        if row == 1
            title(ax, comps{col}, 'FontSize',16,'FontWeight','bold');
        end

        if col == 1
            if(ax.YAxisLocation == "right")
                ax.YTickLabel = [];
                ax.YAxisLocation = "left";
            end
            ylabel(ax, GPSNameList{iSite}, ...
               'Rotation',90, 'FontSize',18, 'FontWeight','bold', ...
               'HorizontalAlignment','center', 'Color', label_color);

        elseif ax.YAxisLocation == "left"
            ax.YTickLabel = [];
        end

        box(ax,'on');
        ax.LineWidth = 1.5;  % or 2 for a heavier border

        ylim(ax, yLims(row,:));
        hold(ax,'off');
    end
end

% —— Append SDH row: two components + legend in the last tile —— 
if hasSDH
    row = nPlot + 1;

    % Columns 1–2: East & North from tiltx/tilty; band from last-comp u_low/u_high
    for col = 1:2
        ax = nexttile((row-1)*3 + col); hold(ax,'on');

        low   = squeeze(u_low(end, col, 1:end-finalindex));
        high  = squeeze(u_high(end, col, 1:end-finalindex));
        ts_row  = ts(:)'; low_row = low(:)'; hi_row = high(:)';
        x_fill = [ts_row, fliplr(ts_row)];
        y_fill = [low_row, fliplr(hi_row)];
        fill(ax, x_fill, y_fill, fitCol, ...
             'EdgeColor','none', 'FaceAlpha',0.6, 'HandleVisibility','off');

        if col == 1
            d = dE; m = mE;
        else
            d = dN; m = mN;
        end
        plot(ax, ts, d, '-', 'Color', dataCol, 'LineWidth',1.2);
        plot(ax, ts, m, '-', 'Color', fitCol,  'LineWidth',1.6);

        xline(ax, 0, 'k:','LineWidth',1);
        grid(ax,'on'); ax.FontSize = 16;

        if col == 1
            ylabel(ax,'SDH (µrad)', 'Rotation',90, 'FontSize',18, 'FontWeight','bold', ...
                   'HorizontalAlignment','center');
        else
            ax.YTickLabel = [];
        end
        % ax.XTickLabel = [];

        % SDH row sits below the original grid; no need to alter above rows' labels
        xlabel(ax,'','FontSize',10);
        box(ax,'on');
        ax.LineWidth = 1.5;  % or 2 for a heavier border

        ylim(ax, yLimsSDH);
        hold(ax,'off');
    end

    % Column 3: legend hosted in the "Up" slot for SDH
    axLeg = nexttile((row-1)*3 + 3); cla(axLeg); hold(axLeg,'on');
    axis(axLeg,'off');
    plot(axLeg, nan, nan, '-', 'Color', dataCol, 'LineWidth',6, ...
         'DisplayName','GPS');
    plot(axLeg, nan, nan, '-', 'Color', fitCol,  'LineWidth',6, ...
         'DisplayName','LSQ fit');
    legend(axLeg,'show', 'Location','north', ...
           'Orientation','vertical', 'Box','off', 'FontSize', 28);
    box(axLeg,'on');
    axLeg.LineWidth = 1.5;  % or 2 for a heavier border
end

% —— Export as high-res PNG ——
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 17 22]); % 10x8 inch canvas
print(gcf,'-dpng','-r200','./PaperFigs/disp_grid.png')



    %% Print out statistics
    
    disp("Net HMM dp: " + (dp(1, 1) - dp(end - finalindex, 1)))
    disp("Net SC dp: " + (dp(1, 2) - dp(end - finalindex, 2)))

    GPSrms = u - permute(usim(:, :, 1:14), [3, 2, 1]);
    GPSrms = GPSrms(:);
    GPSrms = rms(GPSrms, 'omitnan');
    disp("GPS RMS Misfit: " + GPSrms)

    % tiltrms = rms(simtiltx - tiltx) + rms(simtilty - tilty);
    % disp("Tilt Unweighted RMS Misfit: " + tiltrms);

    clear GPSrms tiltrms
    
    %% Load in geometry
    optimizedM = get_full_m(taiyi_parameters, optParams, true, "gps");
    mHMM = optimizedM(1:8);
    mSC = optimizedM(9:end);
    

    %% Quiver Plot
    % [gHMM, ~, ~, ~] = spheroid(mHMM, [x; y; z(1:length(x))], 0.25, 3.08*10^9);
    % [gSC, ~, ~, ~] = spheroid(mSC, [x; y; z(1:length(x))], 0.25, 3.08*10^9);

    [gHMM, gSC] = creategreens(mHMM, mSC);
    [gTiltHMM, gTiltSC] = createtiltgreens(mHMM, mSC, 0, false);

    % x = GPS_llh(1,:);
    % y = GPS_llh(2,:);
    
    pred_color = "#A2142F";
    obs_color = "#026acc";

    tiltscale = 2e-3;
    u1d(end + 1, 1) = tiltreduced(1) .* tiltscale;
    u1d(end, 2) = tiltreduced(2) .* tiltscale;
    u1d(end, 3) = zeros(length(tiltreduced(1)), 1);
    x(end+1) = xtilt;
    y(end+1) = ytilt;
    nanstat(end + 1) = false;
   
    gtot = gHMM + gSC;
    gtot(:, end + 1) = [(gTiltHMM + gTiltSC), 0] .* tiltscale;

    % Set up optimization result vectors from LSQ solution
    u1d_LSQ = usim(end-finalindex, :, :);
    u1d_LSQ = squeeze(u1d_LSQ(1, :, :))';
    u1d_LSQ(end, 1) = u1d_LSQ(end, 1) .* tiltscale;
    u1d_LSQ(end, 2) = u1d_LSQ(end, 2) .* tiltscale;

    % Delete tilt (temporary)
    % u1d_LSQ(end,:) = [];
    
    hold off;
    figquiver = figure(6);
    clf;
    hold on;
    nanstat_quiver = nanstat;
    nanstat_quiver(end) = true;
    gps_ind = true(size(u1d, 1), 1);
    gps_ind(end) = false;
    sim_ind = true(length(x), 1); sim_ind(end) = false;

    % Create vectors for horizontal displacements GPS stations
    realquiver = quiver3(x(~nanstat_quiver)', y(~nanstat_quiver)', zeros(size(x(~nanstat_quiver)))', u1d(gps_ind, 1) * radscale, u1d(gps_ind, 2) * radscale, ...
        u1d(gps_ind, 3) * 0, 'AutoScale', 'off', 'LineWidth',2, 'MaxHeadSize', 0.1, 'Color', obs_color, 'DisplayName', 'Data');
    simquiver = quiver3(x(~nanstat_quiver)', y(~nanstat_quiver)', zeros(size(x(~nanstat_quiver)))', u1d_LSQ(~nanstat_quiver, 1) * radscale, u1d_LSQ(~nanstat_quiver, 2) * radscale, ...
        u1d_LSQ(~nanstat_quiver, 3) * 0, 'AutoScale', 'off', 'LineWidth',2, 'MaxHeadSize', 0.1, 'Color', pred_color, 'DisplayName','Optimization Result');

    % Create tilt quiver plot
    quiver3(x(end)', y(end)', zeros(size(x(end)))', u1d(end, 1) * radscale, u1d(end, 2) * radscale, ...
        u1d(end, 3) * 0, 'AutoScale', 'off', 'LineWidth',5, 'MaxHeadSize', 0.3, 'Color', obs_color, 'HandleVisibility','off');
    quiver3(x(end)', y(end)', zeros(size(x(end)))', u1d_LSQ(end, 1) * radscale, u1d_LSQ(end, 2) * radscale, ...
        u1d_LSQ(end, 3) * 0, 'AutoScale', 'off', 'LineWidth',5, 'MaxHeadSize', 0.3, 'Color', pred_color, 'HandleVisibility','off');
    
    % GNSS stations: dots to show location
    plot3(x(~nanstat_quiver), y(~nanstat_quiver), zeros(size(y(~nanstat_quiver))), '.', ...
          'Color','k','MarkerSize',20,'HandleVisibility','off');
    
    % Tilt-meter: draw a filled triangle for location
    plot3(x(end), y(end), 0, ...
          '^','MarkerSize',15,'MarkerFaceColor','k',...
          'MarkerEdgeColor','k','HandleVisibility','off');

    % ONLY PLOT TILT
    % realquiver = quiver3(x(end)', y(end)', zeros(size(x(end)))', u1d(end, 1) * radscale, u1d(end, 2) * radscale, ...
    %     u1d(end, 3) * 0, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#4DBEEE', 'DisplayName', 'Data');
    % simquiver = quiver3(x(end)', y(end)', zeros(size(x(end)))', u1d_LSQ(end, 1) * radscale, u1d_LSQ(end, 2) * radscale, ...
    %     u1d_LSQ(end, 3) * 0, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#A2142F', 'DisplayName','Optimization Result');
    

    % resquiver = quiver3(x(~nanstat)', y(~nanstat)', zeros(size(x(~nanstat)))', (u1d(:, 1)-u1d_LSQ(~nanstat, 1)) * radscale, (u1d(:, 2)-u1d_LSQ(~nanstat, 2)) * radscale, ...
    %     (u1d(:, 3) - u1d_LSQ(~nanstat, 3)), 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#4DBEEE', 'DisplayName', 'Residual');
    % Plot vert. displacement circles:
    theta = linspace(0, 2*pi, 50);  
    verticalScale = 0.1e4; 
    u1d_ind = 1;
    for k = 1:length(u1d_LSQ)
        if(~nanstat(k))
            line_style = '--';
            if(u1d(u1d_ind,3) < 0); line_style = '-'; end

            r = abs(u1d(u1d_ind,3)) * verticalScale;
            x_data = x(k) + r*cos(theta);
            y_data = y(k) + r*sin(theta);

            plot(x_data, y_data, line_style, 'Color', obs_color, 'LineWidth', 3, 'HandleVisibility','off');
            u1d_ind = u1d_ind + 1;
        
            % Take this out of if statement if we want to plot pred. for
            % simulated data
            line_style = '--';
            if(u1d_LSQ(k,3) < 0); line_style = '-'; end
    
            r_lsq = abs(u1d_LSQ(k,3)) * verticalScale;
    
            x_sim = x(k) + r_lsq*cos(theta);
            y_sim = y(k) + r_lsq*sin(theta);
    
            plot(x_sim, y_sim, line_style, 'Color', pred_color, 'LineWidth', 2.5, 'HandleVisibility','off');
            axis equal;
        end
    end

    % Add reference circle for scale
    % Define reference location and scale
    x_ref = 6000;
    y_ref = 6000;
    r_ref = 0.3 * verticalScale;   % Reference vertical displacement magnitude (e.g., 1 unit)
    x_ref_circle = x_ref + r_ref*cos(theta);
    y_ref_circle = y_ref + r_ref*sin(theta);
    x_ref_tilt = 6e3;
    y_ref_tilt = 5e3;
    
    % Plot the reference circle
    plot(x_ref_circle, y_ref_circle, 'k', 'LineWidth', 3, 'HandleVisibility','off');
    % text(x_ref - 250, y_ref - r_ref - 200, "30 cm", 'Color', 'k', 'FontSize', 20, 'HandleVisibility','off');
    quiver(x_ref, y_ref, -radscale * 0.5, 0, 'LineWidth',4, 'MaxHeadSize', 0.6, 'Color', 'k', 'HandleVisibility','off'); % Vector scale reference
    % text(x_ref - radscale*0.4, y_ref -radscale*0.1, "50 cm", 'Color', 'k', 'FontSize', 20);

    % Add reference vector for SDH tiltmeter
    quiver(x_ref_tilt, y_ref_tilt, -tiltscale*radscale*100, 0, 'LineWidth',4, 'MaxHeadSize', 0.6, 'Color', 'k', 'HandleVisibility','off')
    % text(x_ref_tilt - radscale*0.15, y_ref_tilt -radscale*0.1, "100 µrad", 'Color', 'k', 'FontSize', 20);

    % Plot ellipsoids
    for i = 1:2
        if(i == 1); m = mHMM; name = "Magma Chamber"; col = "#ffa07d"; end
        if(i == 2); m = mSC; name = "SC"; col = "#ffa07d"; end %#ff7d7d

        % 1.  Ellipsoid parameters
        x_c = m(5);   y_c = m(6);   z_c = m(7);
        a   = m(2);                       %  x-semi-axis
        b   = m(1);                       %  y-semi-axis
        c   = m(2);                       %  z-semi-axis  (update to your data)
        
        % 2.  Build a sphere mesh, then rotate so the +z axis is the 
        %     ellipsoid’s longest semi-axis.  That gives us the right
        %     “north-pole” for the grid
        n = 20;                                 % mesh resolution
        [phi,theta] = meshgrid(linspace(0,2*pi,n), linspace(0,pi,n));
        xs =  sin(theta) .* cos(phi);           % unit sphere
        ys =  sin(theta) .* sin(phi);
        zs =  cos(theta);
        
        % Which axis is longest?
        [~,imax] = max([a, b, c]); % 1=x, 2=y, 3=z
        if(i == 2); imax = 2; end
        switch imax  % rotate sphere → +z
            case 1  % a is longest  →  rotate −90° about y
                R_align = [ 0  0  1;
                             0  1  0;
                            -1  0  0];
            case 2 % b is longest  →  rotate +90° about x
                R_align = [ 1  0  0;
                             0  0  1;
                             0 -1  0];
            otherwise  % c already longest → no rotation
                R_align = eye(3);
        end
        Sph = R_align * [xs(:)'; ys(:)'; zs(:)'];      % 3 × (N²)
        
        xs = reshape(Sph(1,:), size(xs));
        ys = reshape(Sph(2,:), size(xs));
        zs = reshape(Sph(3,:), size(xs));
        
        % 3.  Scale to the requested semi-axes
        x_ell = a * xs;      y_ell = b * ys;      z_ell = c * zs;
        
        % 4.  Apply strike / dip rotations and translation
        strike = deg2rad(m(4));                   % around z
        dip    = deg2rad(m(3));                   % around x
        Rz = [ cos(strike)  sin(strike)  0;
              -sin(strike)  cos(strike)  0;
                          0           0  1];
        Rx = [ 1           0            0;
               0  cos(dip) -sin(dip);
               0  sin(dip)  cos(dip)];
        Rot = Rx * Rz;
        
        pts = Rot * [x_ell(:)'; y_ell(:)'; z_ell(:)'];
        x_rot = reshape(pts(1,:), size(x_ell)) + x_c;
        y_rot = reshape(pts(2,:), size(y_ell)) + y_c;
        z_rot = reshape(pts(3,:), size(z_ell)) + z_c;
        
        % ------------------------------------------------------------------
        % 5.  Plot with a visible grid
        % ------------------------------------------------------------------
        if(i==1); surf(x_rot, y_rot, z_rot, 'FaceColor', col, 'FaceAlpha', 0.8, 'EdgeColor', 'k', 'EdgeAlpha', 0.5, 'DisplayName', name);
        else; surf(x_rot, y_rot, z_rot, 'FaceColor', col, 'FaceAlpha', 0.8, 'EdgeColor', 'k', 'EdgeAlpha', 0.5, "HandleVisibility", 'off'); end
        hold on

    end

    % Add lighting effects
    light('Position', [0, 0, 200], 'Style', 'infinite'); % Add light source
    % lighting gouraud; % Smooth shading
    % material shiny; % Enhance reflectivity for a polished appearance
    axis equal
    lighting gouraud
    camlight headlight

    addScaleBar(gca, 2000, 8, 150, ...
            'Color1',[.3 .3 .3], 'Color2',[.8 .8 .8], 'FaceAlpha',1);

    % xlabel('Easting (m)', "FontSize", 18);
    % ylabel('Northing (m)', "FontSize", 18);
    % zlabel('Vertical (m)',  "FontSize", 18);
    % title("Displacements at various stations (y oriented N, x oriented E)", "FontSize", 18)
    xlim([-7.1e3, 7.1e3]);
    ylim([-7.1e3, 7.1e3]);
    zlim([-7.1e3, 7.1e3]);
    GPSNameList(end + 1) = "SDH";
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    
    cxy = llh2local(coast_new', [-155.2784, 19.4073]);
    cxy = cxy * 1000;
    plot(cxy(1, :)', cxy(2, :)', 'k-', 'HandleVisibility','off');
    
    % text((x(~nanstat)-300)', (y(~nanstat)+300)', (zeros(size(x(~nanstat))) + 400)', GPSNameList(~nanstat), 'FontSize', 16, 'HandleVisibility','off');
    hold off;
    legend("FontSize", 28, "Location", "southwest");

    if(saveFigs); exportgraphics(figquiver, './PaperFigs/quiver_plot.png', 'Resolution', 500); end


    %% Plot just the reservoir geometry:
        %% Quiver Plot
    
    hold off;
    figquiver = figure(7);
    clf;
    hold on;
    
    % GNSS stations: dots to show location
    plot3(x.*1e-3, y.*1e-3, zeros(size(y)), '.', ...
          'Color','#bf0000','MarkerSize',20,'HandleVisibility','off');
    
    % Tilt-meter: draw a filled triangle for location
    plot3(x(end).*1e-3, y(end).*1e-3, 0, ...
          '^','MarkerSize',15,'MarkerFaceColor','#bf0000',...
          'MarkerEdgeColor','#bf0000','HandleVisibility','off');

   
    % Plot ellipsoids
    for i = 1:2
        if(i == 1); m = mHMM; name = "Magma Chamber"; col = "#ffa07d"; end
        if(i == 2); m = mSC; name = "SC"; col = "#ffa07d"; end %#ff7d7d

        % 1.  Ellipsoid parameters
        x_c = m(5);   y_c = m(6);   z_c = m(7);
        a   = m(2);                       %  x-semi-axis
        b   = m(1);                       %  y-semi-axis
        c   = m(2);                       %  z-semi-axis  (update to your data)
        
        % 2.  Build a sphere mesh, then rotate so the +z axis is the 
        %     ellipsoid’s longest semi-axis.  That gives us the right
        %     “north-pole” for the grid
        n = 20;                                 % mesh resolution
        [phi,theta] = meshgrid(linspace(0,2*pi,n), linspace(0,pi,n));
        xs =  sin(theta) .* cos(phi);           % unit sphere
        ys =  sin(theta) .* sin(phi);
        zs =  cos(theta);
        
        % Which axis is longest?
        [~,imax] = max([a, b, c]); % 1=x, 2=y, 3=z
        if(i == 2); imax = 2; end
        switch imax  % rotate sphere → +z
            case 1  % a is longest  →  rotate −90° about y
                R_align = [ 0  0  1;
                             0  1  0;
                            -1  0  0];
            case 2 % b is longest  →  rotate +90° about x
                R_align = [ 1  0  0;
                             0  0  1;
                             0 -1  0];
            otherwise  % c already longest → no rotation
                R_align = eye(3);
        end
        Sph = R_align * [xs(:)'; ys(:)'; zs(:)'];      % 3 × (N²)
        
        xs = reshape(Sph(1,:), size(xs));
        ys = reshape(Sph(2,:), size(xs));
        zs = reshape(Sph(3,:), size(xs));
        
        % 3.  Scale to the requested semi-axes
        x_ell = a * xs;      y_ell = b * ys;      z_ell = c * zs;
        
        % 4.  Apply strike / dip rotations and translation
        strike = deg2rad(m(4));                   % around z
        dip    = deg2rad(m(3));                   % around x
        Rz = [ cos(strike)  sin(strike)  0;
              -sin(strike)  cos(strike)  0;
                          0           0  1];
        Rx = [ 1           0            0;
               0  cos(dip) -sin(dip);
               0  sin(dip)  cos(dip)];
        Rot = Rx * Rz;
        
        pts = Rot * [x_ell(:)'; y_ell(:)'; z_ell(:)'];
        x_rot = reshape(pts(1,:), size(x_ell)) + x_c;
        y_rot = reshape(pts(2,:), size(y_ell)) + y_c;
        z_rot = reshape(pts(3,:), size(z_ell)) + z_c;
        
        % ------------------------------------------------------------------
        % 5.  Plot with a visible grid
        % ------------------------------------------------------------------
        if(i==1); surf(x_rot.*1e-3, y_rot.*1e-3, z_rot.*1e-3, 'FaceColor', col, 'FaceAlpha', 0.8, 'EdgeColor', 'k', 'EdgeAlpha', 0.5, 'DisplayName', name);
        else; surf(x_rot.*1e-3, y_rot.*1e-3, z_rot.*1e-3, 'FaceColor', col, 'FaceAlpha', 0.8, 'EdgeColor', 'k', 'EdgeAlpha', 0.5, "HandleVisibility", 'off'); end
        hold on

    end

    % Add lighting effects
    light('Position', [0, 0, 200], 'Style', 'infinite'); % Add light source
    axis equal
    lighting gouraud
    camlight headlight

    xlabel('Easting (km)', "FontSize", 20);
    ylabel('Northing (km)', "FontSize", 20);
    zlabel('Vertical (km)',  "FontSize", 20);
    % title("Displacements at various stations (y oriented N, x oriented E)", "FontSize", 18)
    xlim([-7.1, 7.1]);
    ylim([-7.1, 7.1]);
    zlim([-7.1, 7.1]);
    % set(gca, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
    
    cxy = llh2local(coast_new', [-155.2784, 19.4073]);
    cxy = cxy * 1000;
    plot(cxy(1, :)'.*1e-3, cxy(2, :)'.*1e-3, 'k-', 'HandleVisibility','off');
    ax = gca;
    ax.FontSize = 18;
    % text((x-300)'.*1e-3, (y+300)'.*1e-3, (zeros(size(x)) + 400)'.*1e-3, GPSNameList, 'FontSize', 20, 'HandleVisibility','off');
    hold off;
    % legend("FontSize", 18, "Location", "southwest");

    if(saveFigs); exportgraphics(figquiver, './PaperFigs/spheroid_plot.png', 'Resolution', 500); end
    %% Plot inSAR
    % cLimits = [-1.0, 1.0];
    % cmap = parula;
    % figinsar = figure(8);
    % 
    % clf;
    % tl2 = tiledlayout(2,2,'Padding','compact', 'TileSpacing','compact');
    % nexttile(tl2);
    % plot_insar(insarx, insary, insaru, block_size, look, x, y, u1d, u1d_LSQ, nanstat, ...
    % radscale, 31, GPSNameList, coast_new,cLimits, cmap, saveFigs);
    % title("Insar data");
    % nexttile(tl2);
    % plot_insar(insarx, insary, insaru_pred, block_size, look, x, y, u1d, u1d_LSQ, nanstat, ...
    % radscale, 31, GPSNameList, coast_new,cLimits, cmap, saveFigs);
    % title("Insar prediction")
    % nexttile(tl2);
    % 
    % cLimits = [-0.5, 0.5];
    % cmap = brewermap(64, 'RdYlBu');
    % plot_insar(insarx, insary, insaru - insaru_pred', block_size, look, x, y, u1d, u1d_LSQ, nanstat, ...
    % radscale, 31, GPSNameList, coast_new, cLimits, cmap, saveFigs);
    % title("Residual");
    % 
    % if(saveFigs); exportgraphics(figinsar, './PaperFigs/insarfit.png', 'Resolution', 500); end
end

function addScaleBar(ax, totalLen, nSeg, width, varargin)
% addScaleBar(ax, totalLen, nSeg, width, Name,Value,…)
%
%   ax        = target axes (gca if omitted)
%   totalLen  = bar length in data units   (m)
%   nSeg      = number of alternating boxes
%   width     = bar height / thickness     (m)
%
% Name-value options:
%   'FaceAlpha'  – transparency (default 0.9)
%   'Color1'     – colour of the first box (default 'k')
%   'Color2'     – colour of the second box (default 'w')
%   'Text'       – label string (default sprintf('%g km',totalLen/1000))

if nargin < 1 || isempty(ax),  ax = gca;  end
p = inputParser;
addParameter(p, 'FaceAlpha', 0.9);
addParameter(p, 'Color1',    'k');     % dark box colour
addParameter(p, 'Color2',    'w');     % light box colour
addParameter(p, 'Text',      sprintf('%.0f km', totalLen/1000));
parse(p, varargin{:});
opt = p.Results;

%% --- choose anchor in *upper-left* corner ------------------------------
xl = ax.XLim;  yl = ax.YLim;  zl = ax.ZLim;

x0 = xl(1) + 0.05*range(xl);               % 5 % from left edge
y0 = yl(2) - 0.05*range(yl) - width;       % 5 % down from top edge
z0 = zl(2) + 0.02*range(zl);               % float *above* top Z

segLen = totalLen/nSeg;                    % each box length

%% --- draw the checkerboard ---------------------------------------------
for k = 0:nSeg-1
    x1 = x0 + k*segLen;
    x2 = x1 + segLen;
    verts = [x1 y0 z0;  x2 y0 z0;  x2 y0+width z0;  x1 y0+width z0];

    patch('Vertices',verts, 'Faces',[1 2 3 4], ...
          'Parent',ax, ...
          'FaceColor',    opt.(sprintf('Color%d',mod(k,2)+1)), ...
          'FaceAlpha',    opt.FaceAlpha, ...
          'EdgeColor',    'k', 'LineWidth',0.5, ...
          'Clipping','off', 'HandleVisibility','off');              % keep it visible on top
end

%% --- add the label ------------------------------------------------------
text(x0 + totalLen/2, y0 + 1.4*width, z0, opt.Text, ...
     'Parent', ax, ...
     'HorizontalAlignment','center', ...
     'FontWeight','bold', ...
     'Clipping','off', 'HandleVisibility','off');                   % label should float too
end
