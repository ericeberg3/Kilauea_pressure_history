function make_disp_plot(t, finalindex, GPSNameList, ux, uy, uz, usim, u_low, u_high, tiltx, tilty, includeInactive)
% PLOTGPSGRID Generates a tiled grid of GPS displacements and Tilt/SDH data.
%
% INPUTS:
%   t             : Time vector
%   finalindex    : Integer index to truncate data from the end
%   GPSNameList   : String array of station names
%   ux, uy, uz    : Observed displacement arrays [nStations x nTime]
%   usim          : Modeled displacement array [nTime x 3 x nStations]
%   u_low, u_high : Confidence interval arrays [nStations x 3 x nTime]
%   tiltx, tilty  : Observed tilt vectors
%   includeInactive : Boolean (true/false) to show UWEV, BYRL, CRIM

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
end