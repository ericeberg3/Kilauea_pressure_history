function plot_insar_new(insarx, insary, insaru, block_size, look, ...
        x, y, u1d, u1d_LSQ, xon, yon, con, grid_res, ...
        GPSNameList, optimizedM, coast_new, cLimits, opacity, saveFigs, dir)
% PLOT_INSAR_NEW — Draw InSAR quadtree on a hill-shaded grey DEM backdrop.
%
% New optional arg #17  demFile : full path to GeoTIFF DEM
%                        (defaults './Data/kilauea_DEM.tif')
% If the file is absent a synthetic hill-shaded relief is generated.

if(dir == "asc") look = look(:,2);
else; look = look(:,1); end

fig = figure(8); clf;

% ----------------------------------------------------------------------
% FILES – edit these two paths
% ----------------------------------------------------------------------
tifFile = "./Data/USGS_13.tif";
pngFile = "./Data/dem_wide.png";

% ----------------------------------------------------------------------
% CONSTANTS
% ----------------------------------------------------------------------
lat0  = 19.4073;          % Kīlauea summit  (origin for ENU)
lon0  = -155.2784;
R_E   = 6371000;          % spherical Earth radius  (m)
sunAz = 315;              % hill-shade azimuth (° clockwise from N)
sunEl = 45;               % hill-shade elevation (° above horizon)

%% ======================================================================
% (1)  LOW-RES GeoTIFF   →   hill-shade backdrop
% ======================================================================
[Zlo,Rlo] = readgeoraster(tifFile,'OutputType','double');
infoLo = georasterinfo(tifFile);
if ~isempty(infoLo.MissingDataIndicator)
    Zlo(Zlo == infoLo.MissingDataIndicator) = NaN;
end

% ─ coordinates ----------------------------------------------------------
if isa(Rlo,'map.rasterref.MapCellsReference')        % PROJECTED GeoTIFF
    cellX = Rlo.CellExtentInWorldX;
    cellY = abs(Rlo.CellExtentInWorldY);

    xv = Rlo.XWorldLimits(1) + cellX*(0:size(Zlo,2)-1) + cellX/2;

    if strcmpi(Rlo.ColumnsStartFrom,'north')
        y0 = Rlo.YWorldLimits(2) - cellY/2;
        yv = y0 - (0:size(Zlo,1)-1)*cellY;      % north→south
        Zlo = flipud(Zlo);  yv = fliplr(yv);    % make yv ascending (north-up)
    else
        yv = Rlo.YWorldLimits(1) + cellY*(0:size(Zlo,1)-1) + cellY/2;
    end

    if ~isempty(Rlo.ProjectedCRS)
        [x0proj,y0proj] = projfwd(Rlo.ProjectedCRS,lat0,lon0);
        xv = xv - x0proj;  yv = yv - y0proj;    % local ENU (m)
        
    end
else                                                 % GEOGRAPHIC GeoTIFF
    cellLon = Rlo.CellExtentInLongitude;
    cellLat = Rlo.CellExtentInLatitude;

    lonv = Rlo.LongitudeLimits(1) + cellLon*(0:size(Zlo,2)-1) + cellLon/2;
    latv = Rlo.LatitudeLimits(2) - cellLat/2 - cellLat*(0:size(Zlo,1)-1);
    Zlo  = flipud(Zlo);  latv = fliplr(latv);   % north-up

    xv = deg2rad(lonv - lon0).*cosd(lat0)*R_E;
    yv = deg2rad(latv - lat0)             *R_E;

    cellX = abs(diff(xv(1:2)));
    cellY = abs(diff(yv(1:2)));
end

% ─ hill-shade -----------------------------------------------------------
az = deg2rad(sunAz);  el = deg2rad(sunEl);
[dzdx,dzdy] = gradient(Zlo,cellX,cellY);
slope  = atan(hypot(dzdx,dzdy));
aspect = atan2(dzdy,-dzdx);
hsLo   = cos(el).*cos(slope) + sin(el).*sin(slope).*cos(az-aspect);
hsLo   = mat2gray(hsLo).^1.3;             % gamma boost
bgImg  = hsLo;                            % start composite with coarse HS


xlow = -1e4; xhigh = 1e4; 
[~, xlow_ind] = min(abs(xv - xlow)); [~, xhigh_ind] = min(abs(xv - xhigh)); 
ylow = -1e4; yhigh = 1e4;
[~, ylow_ind] = min(abs(yv - ylow)); [~, yhigh_ind] = min(abs(yv - yhigh)); 

% Crop coordinates to more manageable bounds
xv_crop = xv(xlow_ind:xhigh_ind); yv_crop = yv(ylow_ind:yhigh_ind);
lat_crop = latv(ylow_ind:yhigh_ind); lon_crop = lonv(xlow_ind:xhigh_ind);

bgImg_crop = bgImg(ylow_ind:yhigh_ind, xlow_ind:xhigh_ind);

background_ax = axes('Parent',fig);
colormap(background_ax, "gray");
imagesc(background_ax, 'XData', lat_crop, 'YData', lon_crop, 'CData', bgImg_crop);
axis(background_ax,'image','off');     % tight, equal axes
hold(background_ax,'on');

%% -----------------------------------------------------------------------
% (2)  InSAR data  – interpolated grid with transparency
% -----------------------------------------------------------------------

% --- 2A · make sure inputs are column vectors ---------------------------
insarx     = insarx(:);
insary     = insary(:);
insaru     = insaru(:);
block_size = block_size(:);          % one entry per quadtree pixel

% --- 2B · interpolate onto the DEM grid ---------------------------------
F = scatteredInterpolant(insarx,insary,insaru,'natural');
[X,Y] = meshgrid(xv_crop,yv_crop);           % same grid that the DEM is on
U_grid = F(X,Y);
cmax = 1.5;
epsPix = 1;

% --- 2C · build a logical mask that follows the quadtree layout ---------
%
% Every InSAR sample is the *centre* of a square whose half-width is
%   ½·block_size(k)·grid_res   (metres).
% Everything inside at least one of those squares is “valid data”.
%
alpha_map         = false(numel(yv_crop),numel(xv_crop));   % start fully transparent
halfSide_m        = 0.5*block_size.*grid_res;     % half-width in metres
extra             = epsPix;                       % small overlap to hide seams

for k = 1:numel(insarx)
    % physical bounds of the k-th quadtree block
    xlo = insarx(k) - halfSide_m(k) - extra;
    xhi = insarx(k) + halfSide_m(k) + extra;
    ylo = insary(k) - halfSide_m(k) - extra;
    yhi = insary(k) + halfSide_m(k) + extra;

    % columns (x) and rows (y) of grid points that lie inside that square
    xMask = xv_crop >= xlo & xv_crop <= xhi;           % 1×Nx
    yMask = yv_crop >= ylo & yv_crop <= yhi;           % 1×Ny

    % logical outer product → Ny×Nx; OR-accumulate into global mask
    alpha_map(yMask, xMask) = true;
end

% make data outside the mask NaN so they do not influence the colour axis
U_grid(~alpha_map) = NaN;

inax = axes('Parent',fig, ...
            'Position', background_ax.Position, ...
            'Color','none');     % <- key: transparent axes
% --- 2D · plot the result ------------------------------------------------
h_insar = imagesc(inax, 'XData', lat_crop, 'YData', lon_crop, 'CData', U_grid);
set(h_insar,'AlphaData',alpha_map .* opacity)


%% -----------------------------------------------------------------------
% (3)  Extras – coastline etc. (unchanged)
% ------------------------------------------------------------------------
cxy = llh2local(coast_new',[-155.2784,19.4073])*1000;
% plot(background_ax, cxy(1,:), cxy(2,:),'k.','HandleVisibility','off');

set(inax, 'YDir','normal', ...
          'DataAspectRatio',[1 1 1])
set([background_ax, inax], 'YDir','normal');

metersPerDegLat = R_E*pi/180;
metersPerDegLon = metersPerDegLat*cosd(lat0);
daspect(background_ax, [metersPerDegLon metersPerDegLat 1]);
daspect(inax,          [metersPerDegLon metersPerDegLat 1]);

%% Plot quiver to show look direction

% Define arrow vector in local coordinates
% For example: arrow 500 m to the East, 500 m to the North
u = look(1) * 4e-2;
v = look(2) * 4e-2;
% For asc use -1.4*u and swap text translation
if(dir == "asc"); x0 = 19.35;
else; x0 = 19.35 - 1.4*u; end
y0 = -155.23;

% Plot HMM and SC coords
hold(inax, 'on');
deg_per_rad = 180/pi;
lat = lat0 + (optimizedM(5) / R_E) * deg_per_rad;
lon = lon0 + (optimizedM(6) / (R_E * cosd(lat0))) * deg_per_rad;
plot(inax, lat, lon, 'square', 'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerSize', 22); 

deg_per_rad = 180/pi;
lat = lat0 + (optimizedM(5+8) / R_E) * deg_per_rad;
lon = lon0 + (optimizedM(6+8) / (R_E * cosd(lat0))) * deg_per_rad;
plot(inax, lat, lon, 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerSize', 20); 

% Arrow shaft
plot(inax, [x0, x0+u], [y0, y0+v], 'k-', 'LineWidth', 6);

% Arrowhead as filled triangle
arrowLength = sqrt(u^2 + v^2);
headLength = 0.3 * arrowLength;   % scale head relative to arrow length

% Arrowhead base angle
theta = atan2(v,u);

% Define local triangle points for head (relative to tip at (0,0))
ah_x = [0, -headLength, -headLength];
ah_y = [0,  +headLength/2, -headLength/2];

% Rotate the triangle to match arrow direction
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
rotated = R * [ah_x; ah_y];

% Translate to arrow tip
tip_x = x0 + u*1.05; 
tip_y = y0 + v*1.05;
patch(inax, rotated(1,:) + tip_x, rotated(2,:) + tip_y, 'k', 'EdgeColor','none');

% Add LOS text at the tip
if(dir == "asc")
    text(inax, tip_x+5.5e-3, tip_y, 'LoS', ...
         'Color', 'k', ...
         'FontSize', 42, ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle');
else
        text(inax, tip_x-5.5e-3, tip_y, 'LoS', ... 
         'Color', 'k', ...
         'FontSize', 42, ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle');
end

%% --- pick a high-resolution colormap (here: 512-level jet) -------------
nLevels  = 512;
cmap     = turbo(nLevels); % slanCM("bjy", nLevels); %diverging_colormap(nLevels);
colormap(inax, flipud(cmap));

% keep your existing colour limits
clim(inax, cLimits);

% Set consistent limits once
[~, lat_upper] = min(abs(7e3 - xv_crop)); [~, lat_lower] = min(abs(abs(-7e3 - xv_crop)));
[~, lon_upper] = min(abs(7e3 - yv_crop)); [~, lon_lower] = min(abs(-7e3 - yv_crop));

xlims = [lat_crop(lat_lower), lat_crop(lat_upper)];
ylims = [lon_crop(lon_lower), lon_crop(lon_upper)];

xlim(background_ax, xlims);
ylim(background_ax, ylims);
xlim(inax, xlims);
ylim(inax, ylims);

% Link axes and lock them 
linkaxes([background_ax, inax],'xy');
if(con)
    % regenerate / update the color-bar
    cb = colorbar(inax);
    % cb.Label.String = 'LOS displacement (m)';
    % cb.Location = 'manual';
    % cb.Position = [0.773 0.665 0.04 0.25];
    % cb.TickDirection = 'in';
    % cb.AxisLocation = 'in';
    % cb.Color = 'k';
    cb.FontSize = 42;
end


pos = background_ax.Position;
inax.Position = pos;
set(inax, 'Layer', 'top')
set(gca, 'FontSize', 42, "LineWidth", 3); % 'YTickLabel', []
if(~xon) set(gca,'XTickLabel', []); end
if(~yon) set(gca,'YTickLabel', []); end

if(saveFigs)
    exportgraphics(fig, "./PaperFigs/insar.png", 'Resolution', 500);
end

end
