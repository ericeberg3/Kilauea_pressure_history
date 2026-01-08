function plotDispGrid(t, finalindex, GPSNameList, u_low, u_high, ux, uy, uz, usim, tiltx, tilty, varargin) 
%% Making grid of displacements + tilt

% —— Parse Inputs for Single Plot Mode ——
targetStation = [];
targetCompIdx = [];
baseFont = 18;
baselw = 1.2;

if nargin > 11
    targetStation = varargin{1};
    rawComp = varargin{2};
    baseFont = 40; 
    baselw = 1.4;
    
    % Convert component string to index (1, 2, or 3)
    if ischar(rawComp) || isstring(rawComp)
        if strcmpi(rawComp, 'East'), targetCompIdx = 1; end
        if strcmpi(rawComp, 'North'), targetCompIdx = 2; end
        if strcmpi(rawComp, 'Up'), targetCompIdx = 3; end
    else
        targetCompIdx = rawComp;
    end
end

isSinglePlot = ~isempty(targetStation);

% —— Toggle inclusion of the three inactive stations ——
includeInactive = true; 

% —— Identify which stations to plot ——
if isSinglePlot
    % SINGLE MODE: Find specific station
    if strcmp(targetStation, "SDH")
        plotIdx = []; % Empty GPS list if target is SDH
        hasSDH = true;
    else
        plotIdx = find(matches(GPSNameList, targetStation));
        hasSDH = false;
    end
    
    % Set column range to just the requested component
    colRange = targetCompIdx; 
    
    % Determine Layout Size
    layoutRows = 1;
    layoutCols = 1;
else
    % GRID MODE: Standard logic
    inactive = ismember(GPSNameList,["UWEV","BYRL","CRIM"]);
    allIdx   = 1:numel(GPSNameList);
    if includeInactive
        plotIdx = allIdx;
    else
        plotIdx = allIdx(~inactive);
    end
    hasSDH = true;
    colRange = 1:3;
    
    % Determine Layout Size
    layoutRows = numel(plotIdx) + hasSDH;
    layoutCols = 3;
end

nPlot = numel(plotIdx);

% —— Shared settings ——
comps      = {'East','North','Up'};
ts         = t(1:end-finalindex);         
dataCol    = [0 0.4470 0.7410];           
fitCol     = [0.8500 0.3250 0.0980];      

% —— Precompute y‐limits (Only for active plotIdx) —— 
yLims = zeros(max(plotIdx), 2); % Sized for max index to avoid indexing errors
if nPlot > 0
    for row = 1:nPlot
        iSite = plotIdx(row);
        lo =  inf; hi = -inf;
        
        % We look at all 3 comps to keep Y-lims consistent for the station, 
        % or just the target comp if in single mode.
        if isSinglePlot, loopCols = colRange; else, loopCols = 1:3; end

        for col = loopCols
            low   = squeeze(u_low(iSite, col, 1:end-finalindex));
            high  = squeeze(u_high(iSite, col, 1:end-finalindex));
            switch col
                case 1, dtemp = ux(iSite,1:end-finalindex); mtemp = squeeze(usim(1:end-finalindex,1,iSite));
                case 2, dtemp = uy(iSite,1:end-finalindex); mtemp = squeeze(usim(1:end-finalindex,2,iSite));
                case 3, dtemp = uz(iSite,1:end-finalindex); mtemp = squeeze(usim(1:end-finalindex,3,iSite));
            end
            lo = min([lo; low(:); dtemp(:); mtemp(:)]);
            hi = max([hi; high(:); dtemp(:); mtemp(:)]);
        end
        pad = (hi - lo) * 0.05; 
        yLims(iSite,:) = [lo - pad, hi + pad];
    end
end

% —— SDH y-limits —— 
if hasSDH
    % ... (Existing SDH limit logic, omitted for brevity, logic remains same)
    % Recalculated locally just for the plot
    iTilt = size(u_low,1);
    lowSDH  = squeeze(u_low(end, 3, 1:end-finalindex));
    highSDH = squeeze(u_high(end, 3, 1:end-finalindex));
    dE = tiltx(1:end-finalindex); dN = tilty(1:end-finalindex);
    if ndims(usim) >= 3 && size(usim,3) >= iTilt
        mE = squeeze(usim(1:end-finalindex,1,iTilt)); mN = squeeze(usim(1:end-finalindex,2,iTilt));
    else
        mE = nan(size(ts)); mN = nan(size(ts));
    end
    
    % Only calc limits for requested component if single mode
    vals = [];
    if ismember(1, colRange), vals = [vals; dE(:); mE(:)]; end
    if ismember(2, colRange), vals = [vals; dN(:); mN(:)]; end
    if ismember(3, colRange), vals = [vals; lowSDH(:); highSDH(:)]; end % Up/Band
    
    lo = min(vals); hi = max(vals);
    pad = (hi - lo) * 0.05;
    yLimsSDH = [lo - pad, hi + pad];
end

% —— Create Grid ——
fig = figure(8); clf;
if isSinglePlot
    set(fig, 'Units','normalized','Position',[0.3 0.3 0.4 0.4]); % Smaller figure for single plot
else
    set(fig, 'Units','normalized','Position',[0.1 0.1 0.6 0.8]);
end
tl = tiledlayout(layoutRows, layoutCols, 'TileSpacing','none','Padding','none');

% —— Plot regular GPS stations —— 
for row = 1:nPlot
    iSite = plotIdx(row);
    
    isInactive = ismember(GPSNameList(iSite), ["UWEV","BYRL","CRIM"]);
    if isInactive, label_color = 'r'; else, label_color = 'k'; end

    % Loop over the computed range (1:3 or just targetCompIdx)
    for col = colRange
        ax = nexttile; % Simplified from ((row-1)*3 + col)
        hold(ax,'on');
        
        % Y-Axis location logic (Keep simplistic for single mode)
        if ~isSinglePlot && mod((row-1)*3 + col, 2) == 0
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
            case 1, d = ux(iSite,1:end-finalindex); m = squeeze(usim(1:end-finalindex,1,iSite));
            case 2, d = uy(iSite,1:end-finalindex); m = squeeze(usim(1:end-finalindex,2,iSite));
            case 3, d = uz(iSite,1:end-finalindex); m = squeeze(usim(1:end-finalindex,3,iSite));
        end
        plot(ax, ts, d, '-', 'Color', dataCol, 'LineWidth',baselw);
        plot(ax, ts, m, '-', 'Color', fitCol,  'LineWidth',baselw + 0.4);
        xline(ax, 0, 'k:','LineWidth',1);
        grid(ax,'on'); ax.FontSize = baseFont-2;
        
        if ~isSinglePlot, ax.XTickLabel = []; end
        if row == 1
            title(ax, comps{col}, 'FontSize',baseFont-2,'FontWeight','bold');
        end

        % Y-Label Logic
        if col == 1 || isSinglePlot
             if(ax.YAxisLocation == "right"), ax.YTickLabel = []; ax.YAxisLocation = "left"; end
             ylabel(ax, GPSNameList{iSite}, 'Rotation',90, 'FontSize',baseFont, 'FontWeight','bold', 'HorizontalAlignment','center', 'Color', label_color);
        elseif ax.YAxisLocation == "left"
            ax.YTickLabel = [];
        end

        box(ax,'on'); ax.LineWidth = 1.5; 
        ylim(ax, yLims(iSite,:));
        hold(ax,'off');
    end
end

% —— Append SDH row —— 
if hasSDH
    % Columns 1–2: East & North from tiltx/tilty
    % If single plot, colRange will be 1 or 2. 
    % Note: SDH doesn't usually have an 'Up' plot (col 3 is Legend in grid mode)
    
    sdhLoop = intersect(colRange, [1 2]); % Only plot 1 or 2 for SDH
    
    for col = sdhLoop
        ax = nexttile; hold(ax,'on');

        low   = squeeze(u_low(end, col, 1:end-finalindex));
        high  = squeeze(u_high(end, col, 1:end-finalindex));
        ts_row  = ts(:)'; low_row = low(:)'; hi_row = high(:)';
        x_fill = [ts_row, fliplr(ts_row)];
        y_fill = [low_row, fliplr(hi_row)];
        fill(ax, x_fill, y_fill, fitCol, ...
             'EdgeColor','none', 'FaceAlpha',0.6, 'HandleVisibility','off');

        if col == 1, d = dE; m = mE; else, d = dN; m = mN; end
        plot(ax, ts, d, '-', 'Color', dataCol, 'LineWidth',baselw);
        plot(ax, ts, m, '-', 'Color', fitCol,  'LineWidth',baselw + 0.4);
        xline(ax, 0, 'k:','LineWidth',1);
        grid(ax,'on'); ax.FontSize = baseFont-2;
        
        if col == 1 || isSinglePlot
            ylabel(ax,'SDH (µrad)', 'Rotation',90, 'FontSize',baseFont, 'FontWeight','bold', 'HorizontalAlignment','center');
        else
            ax.YTickLabel = [];
        end
        
        xlabel(ax,'','FontSize',baseFont-8);
        box(ax,'on'); ax.LineWidth = 1.5;
        ylim(ax, yLimsSDH);
        hold(ax,'off');
    end

    % Legend (Only show if in grid mode, or strictly requested? 
    % usually single plots don't need a dedicated legend tile)
    if ~isSinglePlot
        axLeg = nexttile; cla(axLeg); hold(axLeg,'on');
        axis(axLeg,'off');
        plot(axLeg, nan, nan, '-', 'Color', dataCol, 'LineWidth',6, 'DisplayName','GPS');
        plot(axLeg, nan, nan, '-', 'Color', fitCol,  'LineWidth',6, 'DisplayName','LSQ fit');
        legend(axLeg,'show', 'Location','north', 'Orientation','vertical', 'Box','off', 'FontSize', baseFont+10);
    end
end

% —— Export ——
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 17 22]); 
end