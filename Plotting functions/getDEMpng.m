fig = figure(8); clf;

% ----------------------------------------------------------------------
% FILES – edit these two paths
% ----------------------------------------------------------------------
tifFile = '/Users/eric/Desktop/Summer 23 Lab/Matlab/ForEric/no_vol_inversion/Data/USGS_13.tif';

% ----------------------------------------------------------------------
% CONSTANTS
% ----------------------------------------------------------------------
lat0  = 19.4073;          % Kīlauea summit  (origin for ENU)
lon0  = -155.2784;
R_E   = 6371000;          % spherical Earth radius  (m)
sunAz = 315;              % hill-shade azimuth (° clockwise from N)
sunEl = 30;               % hill-shade elevation (° above horizon)
zExagg = 3;
gammaAdj = 1.5;

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

% Apply Vertical Exaggeration
[dzdx,dzdy] = gradient(Zlo * zExagg, cellX, cellY); 
slope  = atan(hypot(dzdx,dzdy));
aspect = atan2(dzdy, -dzdx);

% Calculate raw hillshade
hsLo   = cos(el).*cos(slope) + sin(el).*sin(slope).*cos(az-aspect);

% Normalize to 0-1
hsLo = (hsLo - min(hsLo(:))) / (max(hsLo(:)) - min(hsLo(:)));

% APPLY BRIGHTNESS & CONTRAST
if exist('imadjust', 'file')
    % arg1: image
    % arg2: input range (clip top/bottom 1% for contrast)
    % arg3: output range (empty [] means 0 to 1)
    % arg4: GAMMA (controls brightness)
    hsLo = imadjust(hsLo, stretchlim(hsLo, [0.01 0.99]), [], gammaAdj); 
else
    % Fallback
    hsLo = hsLo.^gammaAdj; 
end

bgImg  = hsLo;

%%
xlow = -1e5; xhigh = 1e5; 
[~, xlow_ind] = min(abs(xv - xlow)); [~, xhigh_ind] = min(abs(xv - xhigh)); 
ylow = -1e5; yhigh = 1e5;
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