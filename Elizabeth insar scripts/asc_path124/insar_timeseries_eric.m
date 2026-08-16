% Eric's Kilauea study site
% Plot for thesis

%% Geographic parameters (stored in dem.rsc)
WIDTH        =  1763      ;    
FILE_LENGTH  =  1403     ;     
X_FIRST      =  -155.35 ;
Y_FIRST      =  19.61   ;
X_STEP       =  0.00027777778;
Y_STEP       =  -0.00027777778;

longitudes = X_FIRST:X_STEP:X_FIRST+X_STEP*(WIDTH-1);
latitudes = Y_FIRST:Y_STEP:Y_FIRST+Y_STEP*(FILE_LENGTH-1);


%% Full-res parameters (stored in elevation.dem.rsc)
fWIDTH       =   8816          ;
fFILE_LENGTH =   2807          ;
fX_FIRST     =   -155.35000000 ;
fY_FIRST     =   19.61000000   ;
fX_STEP      =   0.55555556E-04;
fY_STEP      =   -0.13888889E-03;
flongitudes = fX_FIRST:fX_STEP:fX_FIRST+fX_STEP*(fWIDTH-1);
flatitudes =  fY_FIRST:fY_STEP:fY_FIRST+fY_STEP*(fFILE_LENGTH-1);

% Output image is e, n, u components of vector from ground to satellite 
fid=fopen("S1B_IW_RAW__0SDV_20180508T042948_20180508T043023_010825_013CA6_363B.lookvector",'r');
looks=fread(fid,inf,'double');
fclose(fid); 
load('hawaii_line_new.mat', 'coast_new', 'new_pit')

lookvector = reshape(looks, 3, fWIDTH, fFILE_LENGTH);

east = squeeze(lookvector(1,:,:)).';
north = squeeze(lookvector(2,:,:)).';
up = squeeze(lookvector(3,:,:)).';

% Example plot: showing East look vector. It doesn't vary a ton in the
% study area. 
figure; imagesc(flongitudes, flatitudes, up); colorbar;
set(gca, 'YDir','normal'); axis image; 
cald_cent = [-155.2784, 19.4073];
[~, iLon] = min(abs(flongitudes - cald_cent(1)));
[~, iLat] = min(abs(flatitudes  - cald_cent(2)));

avglook_asc = squeeze( lookvector(:, iLon, iLat) );  
save look_asc.mat avglook_asc;
% avglook_asc = [mean(east(:)); mean(north(:)); mean(up(:))];

%% Get data and find residuals of GPS and Insar
% load insar_pred_griddata.mat;
daily_GPS = readtable('daily_gps_mult_stations.csv');
load('GPS_seism_locations.mat', 'GPSNameList', 'GPS_llh');
nstat = length(GPSNameList); % Number of GPS stations

% station x and y positions
xy = llh2local(GPS_llh, [-155.2784, 19.4073]);
x = xy(1,:) * 1000;
y = xy(2,:) * 1000;
z = zeros(1, nstat - 1);

 %% Realistically, you probably won't need to plot look vectors. You can figure 
 % out your own accuracy, but this code calculates the mean and standard
 % deviation of the ENU. I think at first-order, you're probably okay to
 % use east_mean, north_mean, up_mean for each direction. 
 % Note these look vectors are also for 5/8/18 acquisition. It may be
 % a little different from the others, though also probably not much! 
 east_mean = mean(east(:)); 
 east_std = std(east(:));
 north_mean = mean(north(:)); 
 north_std = std(north(:));
 up_mean = mean(up(:)); 
 up_std = std(up(:));


%% 
% get some parameters
load 'parameters'
nr=parameters(1); % size of image in "range" direction
naz=parameters(2); % size of image in "azimuth" direction
nslcs=parameters(3); % number of SAR acquisitions
nigrams=parameters(4); % number of interferograms

jd=load('jdlist');
t_date=datetime(jd,'convertfrom','juliandate','Format','ddMMMyy','TimeZone','UTC'); % Need to set accurate time of Sentinel overpass: 6am or 6pm on the appropriate days!

intlist = importdata('intlist_sequential');
nigrams = length(intlist);

wrapped_int = zeros(naz,nr,nigrams);
for k = 1:nigrams
    fprintf(strcat('Loading phases_', num2str(k),'...\n'));
    fid=fopen(intlist{k},'r');
    phase=fread(fid,[inf],'float32');
    phasecpx = phase(1:2:end) + 1i*phase(2:2:end);
    clear phase; % clears up space
    wrapped_int(:,:,k)=reshape(phasecpx,nr,naz).'; % Wrapped interferogram
    fclose(fid);
end
%% Coherence load
cc_array = zeros(naz,nr,nigrams);
for k = 1:nigrams
    fprintf(strcat('Loading phases_', num2str(k),'...\n'));
    fid=fopen(strcat(intlist{k}(1:17),'.cc'),'r');
    phase_cc=fread(fid,[nr*2 naz],'float32');
    % clear phase_cc; % clears up space
    cc_array(:,:,k)=phase_cc(nr+1:end,:).'; % Wrapped interferogram
    fclose(fid);
end

%% Loads amplitude to make nice underlying image based on radar brightness


% acc=dat(1:nr,:).';
acc = abs(wrapped_int(:,:,1));
meanacc=mean(acc(:));
stdacc=std(acc(:));
acc(acc>meanacc+1*stdacc)=meanacc+1*stdacc;
acc(acc<meanacc-1*stdacc)=meanacc-1*stdacc;
acc=acc-(max(0,meanacc-3*stdacc));
acc=acc/max(acc(:));
figure;
imagesc(acc.^0.2); colormap gray;
axis image;

fprintf('Interferometric Data loaded.\n');


% Howard colormap
addpath('../../../matlabcodeforfigures/')

%% Mask for acc 
% Mask 
[X,Y] = meshgrid(1:size(acc,1),1:size(acc,2));
mask = 1020*Y - 209*X < 586291; 
% figure; imagesc(disp_amp(:,:,end) + disp_amp(:,:,end).*mask.');
% clim([0 1])
acc = acc + acc .* mask.';


%% Video of interferogram phase
% Test of coherence mask
% test_mask = min(cc_array,[],3) > 0.2; 
% test_mask = mean(cc_array,3)>0.4;
% 
% doPhaseVideo = 1; 
% if doPhaseVideo == 1
% v = VideoWriter('kilauea_timeseries_masktest.mp4','MPEG-4');
% v.FrameRate = 6;
% open(v);
% figure;
% % set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf, 'Position',[160 440 664 584]);
% for k = 1:nigrams
%     rgb = dismph(acc.*test_mask, angle(wrapped_int(:,:,k)), pi,naz, nr);
%     imagesc(longitudes, latitudes, rgb);
%     title(strcat("Phase Difference, ",string(t_date(k))," to " ,string(t_date(k+1))));
%     colorbar;% colormap(howard_colormap); 
%     clim([-pi pi]);
%     set(gca,'YDir','normal')
%     axis image; 
%     xlim([-155.35 -155.1]); ylim([19.22, 19.5]);
% 
%     F = getframe(gcf);
%     writeVideo(v,F);
%     % pause(0.2)
% end
% close(v);
% end

%% Load displacement
fprintf(strcat('Loading Displacement...\n'));
fid=fopen("displacement",'r');
dat=fread(fid,[inf],'float32');
fclose(fid);
disp=reshape(dat,nr*2,naz,nigrams); % Displacement
disp_amp = disp(1:nr,:,:);
displacement = disp(nr+1:end,:,:);
%%
lambda = 5.56e-2;
displacement = displacement/(4*pi)*lambda; % In meters
cum_disp_asc = squeeze(displacement(:,:,end)).';
avg_cc_asc = mean(cc_array, 3);
% max_cc_asc = max(cc_array, 3);
min_cc_asc = min(cc_array, [], 3);
bot10_cc_asc = prctile(cc_array, 10, 3);
bot30_cc_asc = prctile(cc_array, 30, 3);


%%
% Mask 
[X,Y] = meshgrid(1:size(disp_amp(:,:,end),1),1:size(disp_amp(:,:,end),2));
mask = 1020*X - 209*Y < 586291; 
% figure; imagesc(disp_amp(:,:,end) + disp_amp(:,:,end).*mask.');
% clim([0 1])
disp_amp = disp_amp + disp_amp .* mask.';

cmax = 0.5*max(max(max(abs(displacement))));
rgb = dishgt(disp_amp(:,:,end).',displacement(:,:,end).',cmax,naz,nr);
figure(12); imagesc(longitudes, latitudes, rgb);
set(gca,'YDir','normal')
colorbar; clim([-cmax cmax]); 
% figure; imagesc(disp_amp(:,:,1)); 
% figure; imagesc(displacement(:,:,1));

%% Quadtree downsample
% [insaru_asc, blocks_asc] = quad_downsample(cum_disp_asc, avg_cc_asc', longitudes, latitudes, coast_new);

% Adjust tol based on data’s noise and variability characteristics
tol = 0.1;

% Get the original dimensions
[nrows, ncols] = size(cum_disp_asc);

% Define your desired maximum block size
B = 256;

% Compute new dimensions that are multiples of B
new_rows = floor(nrows/B)*B;
new_cols = floor(ncols/B)*B;

% Crop the image
unw_crop = cum_disp_asc(1:new_rows, 1:new_cols);
cc_crop = avg_cc_asc(1:new_rows, 1:new_cols);

cc_tests = zeros(3, size(cc_crop, 1), size(cc_crop, 2));
cc_tests(1, :, :) = cc_crop;% mean cc
cc_tests(2, :, :) = min_cc_asc(1:new_rows, 1:new_cols);
cc_tests(3, :, :) = bot30_cc_asc(1:new_rows, 1:new_cols);


% Perform quadtree decomposition
S = qtdecomp(unw_crop, tol, [1, B]);

% Initialize the resampled image with the original data
resampled = cum_disp_asc;

% Get the list of unique block sizes present in the quadtree decomposition
blockSizes = full(unique(S(S > 0)));
% Process blocks from largest to smallest (you can change the order if needed)
blockSizes = sort(blockSizes, 'descend');

% Loop over each block size
for bs = blockSizes'
    % Extract blocks of the current size using qtgetblk
    blocks = qtgetblk(cum_disp_asc, S, bs);
    if ~isempty(blocks)
        % Find indices (upper left corners) of blocks with size bs in S
        [r, c] = find(S == bs);
        for k = 1:length(r)
            % Define the block indices in cum_disp_asc
            block = cum_disp_asc(r(k):r(k)+bs-1, c(k):c(k)+bs-1);
            % Compute the mean of the block
            blockMean = mean(block(:));
            % Replace the block with its mean value in the resampled image
            resampled(r(k):r(k)+bs-1, c(k):c(k)+bs-1) = blockMean;
        end
    end
end

% Display the original and resampled displacement images
figure(7);
clf;
subplot(1,2,1);
imagesc(longitudes, latitudes, cum_disp_asc);
hold on;
% plot(coast_new(:, 1)', coast_new(:, 2)', 'k.', 'HandleVisibility','off');
colorbar; 
title('Original Displacement (meters)');
set(gca,'YDir','normal');

subplot(1,2,2);
imagesc(longitudes, latitudes, resampled);
hold on;
% plot(coast_new(:, 1)', coast_new(:, 2)', 'k.', 'HandleVisibility','off');
colorbar; 
title('Quadtree Resampled Displacement');
set(gca,'YDir','normal');

%% Filter out bad data:
% Assume these have already been defined and cropped to compatible sizes:
% unw_crop: the cropped displacement data (e.g., 1403x1763 cropped to multiples of maxBlock)
% cc_crop: the cropped coherence data
% S: the quadtree structure from qtdecomp on unw_crop

% Filter out values within caldera
nan_threshold_coh = 0.10;
filtered_resampled = downsample_coh(0.3, nan_threshold_coh, min_cc_asc, resampled, S);
cald_bnds = [[-155.277, 19.43]; [-155.257, 19.41]];

[~, left_ind] = min(abs(longitudes - cald_bnds(1,1)));
[~, right_ind] = min(abs(longitudes - cald_bnds(2,1)));
[~, up_ind] = min(abs(latitudes - cald_bnds(1,2)));
[~, down_ind] = min(abs(latitudes - cald_bnds(2,2)));

filtered_resampled(up_ind:down_ind, left_ind:right_ind) = nan;

function filtered_resampled = downsample_coh(cc_threshold, nan_threshold, cc_crop, resampled, S)
% Define your coherence threshold
% cc_threshold = 0.38;

% Create a copy of your resampled displacement data to apply the coherence filter
filtered_resampled = resampled;

% Get the list of unique block sizes present in the quadtree decomposition (largest to smallest)
blockSizes = full(unique(S(S > 0)));
blockSizes = sort(blockSizes, 'descend');

% Loop through each block size
for bs = blockSizes'
    % Find the indices (upper-left corners) of blocks with size bs
    [r, c] = find(S == bs);
    for k = 1:length(r)
        % Define indices for the current block in both displacement and coherence arrays
        rowIdx = r(k):(r(k)+bs-1);
        colIdx = c(k):(c(k)+bs-1);
        
        % Extract the coherence block from cc_crop
        coherence_block = cc_crop(rowIdx, colIdx);
        % Compute the mean coherence for the block (or you can use another metric)
        % mean_coh = mean(coherence_block(:));
        
        % % Mean coh metric
        % if mean_coh < cc_threshold
        %     filtered_resampled(rowIdx, colIdx) = NaN;
        % end

        % Percent of nan pixels in block metric (Elizabeth suggested)
        coherence_block_filtered = coherence_block < cc_threshold;
        if(sum(coherence_block_filtered) > (nan_threshold * length(coherence_block(:))) )
            filtered_resampled(rowIdx, colIdx) = NaN;
        end
    end
end
end

% Display the filtered displacement data alongside the original resampled data
figure(8);
clf;
subplot(1,2,1);
imagesc(longitudes, latitudes, -resampled);
hold on;
% plot(coast_new(:, 1)', coast_new(:, 2)', 'k.', 'HandleVisibility','off');
colorbar;
title('Quadtree Resampled Displacement');
scatter(GPS_llh(1,:), GPS_llh(2,:), 100, cumDisp_looked, 's', ...
        'filled', 'MarkerEdgeColor','k');
quiver(-155.22, 19.47, avglook_asc(1), avglook_asc(2), 1e-1, "filled", "LineWidth", 4)
text(GPS_llh(1,:) + 0.004, GPS_llh(2,:), GPSNameList, "FontSize", 15)
xlim([-155.35, -155.2]);
ylim([19.35, 19.5]);
set(gca,'YDir','normal');

subplot(1,2,2);
imagesc(longitudes, latitudes, -filtered_resampled);
hold on;
colorbar;
plot(GPS_llh(1,:), GPS_llh(2,:), 'r.', "MarkerSize", 20)
plot(coast_new(:, 1)', coast_new(:, 2)', 'k.', 'HandleVisibility','off');
text(GPS_llh(1,:) + 0.004, GPS_llh(2,:), GPSNameList, "FontSize", 15, "Color", "r")
title('Filtered by Coherence');
xlim([-155.35, -155.2]);
ylim([19.35, 19.5]);
set(gca,'YDir','normal');

%% Set up HMM and SC geometries
% taiyi_parameters = [1600.79, 914.47, 90, 0, 70.3, 183, -2.18e3, -3e6, ... 
%     277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, -10e6];
% taiyi_parameters = [1600.79000000000	914.470000000000	90	0	70.3000000000000	183	-2180	-7578363.95906990 ...
%     453.334521752308	1147.40385217986	67.9911620543601	136	1617.21062995580	-2479.30350291513	-3630	-9224358.48561988];
% 
% 
% XY_grid_m = zeros(size(resampled));
% insaru_pred_HMM = zeros(size(resampled));
% insaru_pred_SC = zeros(size(resampled));
% 
% for i = 1:size(XY_grid_m, 1)
%     parfor j = 1:size(XY_grid_m, 2)
%         xy_temp = llh2local([longitudes(j); latitudes(i)], [-155.2784, 19.4073])'.*1e3;
%         insaru_pred_HMM(i, j) = spheroid(taiyi_parameters(1:8), [xy_temp(1); xy_temp(2); 0], 0.25, 3.08*10^9)' * avglook_asc;
%         insaru_pred_SC(i, j) = spheroid(taiyi_parameters(9:end), [xy_temp(1); xy_temp(2); 0], 0.25, 3.08*10^9)' * avglook_asc;
%     end
% end
% 
% save insar_pred_griddata.mat insaru_pred_SC insaru_pred_HMM

% COMPARING DATES BETWEEN May 08, 2018-August 06, 2018
%%  Assemble daily GPS displacements for all stations, 08-May-2018 → 06-Aug-2018

nStat   = numel(GPSNameList);
cumDisp = nan(nStat,3);  % rows → stations, cols → [E N U]
comps   = ["east","north","up"];
% Wrap in a table for readability
cumDispTable = array2table(cumDisp, ...
           'VariableNames', comps, ...
           'RowNames',      GPSNameList);

cumDisp_looked = cumDisp * avglook_asc;

% 1. Parse the date column --------------------------------------------------
gpsDates = datetime(daily_GPS.Date, ...
                    'InputFormat',"MM/dd/yyyy HH:mm", ...   % works for '05/08/2018 02:00'
                    'TimeZone',"UTC");

% If the file encodes the year as '0018' instead of '2018',
% bump everything with year < 1900 forward by 2000 yr
badYear = year(gpsDates) < 1900;
gpsDates(badYear) = gpsDates(badYear) + calyears(2000);

% 2. Keep only the interval of interest ------------------------------------
t0 = datetime(2018,5,8,'TimeZone','UTC');
t1 = datetime(2018,8,6,'TimeZone','UTC');

inWindow  = gpsDates >= t0 & gpsDates <= t1;
gpsDates  = gpsDates(inWindow);           % dates you keep
nDates    = nnz(inWindow);

t0 = datetime(2018,5,8,'TimeZone','UTC');
t1 = datetime(2018,8,6,'TimeZone','UTC');

% Logical index of rows within the window
inWindow = gpsDates >= t0 & gpsDates <= t1;

% These are the rows we’ll use as “start” and “end”
startRow = find(inWindow,1,'first');
endRow   = find(inWindow,1,'last');

for s = 1:nStat
    base = GPSNameList{s};
    if isstrprop(base(1),'digit') % add “x” if the name begins with a digit
        base = "x" + base;
    end

    for c = 1:3
        varName = base + "_" + comps(c);

        if ismember(varName,daily_GPS.Properties.VariableNames)
            % First non-NaN inside the window
            colData = daily_GPS{inWindow, varName};
            idx_firstval = find(~isnan(colData),1,'first');
            idx_lastval = find(~isnan(colData),1,'last');
            firstVal = colData(idx_firstval);
            lastVal  = colData(idx_lastval);
            days_missed =  ((idx_firstval - 1)+ (length(colData) - idx_lastval));
            fprintf("Days missed at " + varName + ": " + days_missed + "\n")
            if ~isempty(firstVal) && ~isempty(lastVal) && days_missed < 10
                cumDisp(s,c) = lastVal - firstVal;
            else
                cumDisp(s, c) = NaN;
            end
        end
    end
end
%%



% 1. Coordinate vectors (coerce to correct orientation) ---------------------
latVec = latitudes(:); % column
lonVec = longitudes(:).';      % row
resampled_interp = resampled;

% 2. Ensure monotonic increase required by griddedInterpolant --------------
if any(diff(latVec) < 0)            % descending latitude → flip
    latVec    = flipud(latVec);
    resampled_interp = flipud(resampled_interp);
end
if any(diff(lonVec) < 0)            % descending longitude → flip
    lonVec    = fliplr(lonVec);
    resampled_interp = fliplr(resampled_interp);
end

% 3. Build interpolant (linear; NaN outside grid)
F = griddedInterpolant({latVec, lonVec}, resampled_interp, ...
                       'linear', 'none');

% 4. Evaluate grid at exact GPS coordinates
gpsLon = GPS_llh(1,:).';       % nStat×1
gpsLat = GPS_llh(2,:).';

insarAtGPS = F(gpsLat, gpsLon);   % nStat×1 vector

% 5. Residuals:  GPS minus InSAR 
residuals = cumDisp_looked + insarAtGPS;  % nStat×1

%% Display residuals
figure(9);
clf;
r = sqrt(x.^2 + y.^2).*1e-3;
plot(r, residuals, '.', 'MarkerSize', 18);
xlabel('Distance from center of caldera (km)');
ylabel('GPS – InSAR displacement (m)');
title('Cumulative-displacement residuals');
set(groot, 'DefaultAxesFontSize', 18)

hold on
for k = 1:numel(GPSNameList)
    % Slight horizontal offset keeps label from sitting on the marker
    text(r(k)+0.02*range(r), ...      % x-position
         residuals(k), ...            % y-position
         GPSNameList{k}, ...          % label
         'FontSize',14, ...
         'VerticalAlignment','bottom', ...
         'HorizontalAlignment','left');
end
hold off


% filtered_resampled_test = downsample_coh(coh_list(i),nan_list(j), min_cc_asc, resampled, S);
% imagesc(longitudes, latitudes, resampled);
% set(gca, 'YDir', 'normal');
% hold on;
% colormap default;
% colorbar
% plot(coast_new(:, 1)', coast_new(:, 2)', 'k.', 'HandleVisibility','off');
% title("Threshold = " + coh_list(i) + ", nan percent acceptable: " + nan_list(j)*100 + "%");
% xlim([-155.33, -155.22]);
% ylim([19.34, 19.48]);

% Plot GPS and InSAR residuals

%% Test a variety of coh threshold / nan threshold
load insar_pred_griddata.mat
thresh = figure(5);
coh_list = [0.2, 0.3, 0.4];
nan_list = [0.01, 0.05, 0.1];
imRef = -resampled;
ccRef = squeeze(cc_tests(2,:,:));
nrows = length(coh_list);
ncols = 3;
grey_index = size(abyss, 1); 

t = tiledlayout(nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'none');
clf;

for i = 1:nrows
    for j = 1:ncols
        % subplot(nrows, ncols, (i-1)*ncols + j);  % Proper subplot indexing 
        ax = nexttile;
        filtered_resampled_test = downsample_coh(coh_list(i), nan_list(j), ccRef, imRef, S); %min_coh_asc

        % Find all NaN values in the filtered data
        % nan_indices = isnan(filtered_resampled_test);
        % 
        % % Replace NaN values with the index of the grey color in our custom colormap
        % % This is a crucial step to map NaNs to the grey color
        % filtered_resampled_test(nan_indices) = grey_index;

        imagesc(longitudes, latitudes, filtered_resampled_test);
        set(gca, 'YDir', 'normal');
        hold on;
        plot(gpsLon, gpsLat, '.r', "MarkerSize", 14);
        % text(gpsLon + 1e-2, gpsLat, GPSNameList, 'color', 'red');
        colormap bone;
        clim([-1.5, 0.5])
        axis square;
        % colorbar;
        plot(coast_new(:, 1)', coast_new(:, 2)', 'w.', 'HandleVisibility','off');
        % title("Threshold = " + coh_list(i) + ", nan percent acceptable: " + nan_list(j)*100 + "%", "FontSize", 13);
        xlim([-155.33, -155.22]);
        ylim([19.34, 19.48]);
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.FontSize = 22;
cb.Label.String = 'LOS Displacement (m)';
% sgtitle("Real deformation with various filter criteria. Min of all coh slices used", "FontSize", 30)
exportgraphics(thresh, './asc_mask_grid.png', 'Resolution', 500);

%% Plot which GPS stations are excluded

latVec = latitudes(:);  lonVec = longitudes(:).';
if any(diff(latVec)<0),  latVec = flipud(latVec);  imRef = flipud(imRef); end
if any(diff(lonVec)<0),  lonVec = fliplr(lonVec);  imRef = fliplr(imRef); end
Fref   = griddedInterpolant({latVec,lonVec}, imRef,'linear','none');
insarRef = Fref(gpsLat, gpsLon);                 % nStat×1  (m)

nrows    = numel(coh_list);
ncols    = numel(nan_list);

% distance of each GPS station from the caldera centre (km)
r = hypot(x,y) * 1e-3;
meanResMat = nan(nrows,ncols);

figure(25); clf;
for i = 1:nrows
    for j = 1:ncols
        % 2. Build the masked image for this coherence / NaN combination
        filt_im = downsample_coh(coh_list(i), nan_list(j), ccRef, ...
                                 resampled, S);

        % 3. Interpolate masked image at GPS locations
        latV = latitudes(:);  lonV = longitudes(:).';  im2int = filt_im;
        if any(diff(latV)<0),  latV = flipud(latV);  im2int = flipud(im2int); end
        if any(diff(lonV)<0),  lonV = fliplr(lonV);  im2int = fliplr(im2int); end
        Fmask  = griddedInterpolant({latV,lonV}, im2int,'linear','none');
        insarM = Fmask(gpsLat, gpsLon);   % NaN if station lies inside mask

        % 4. GPS – InSAR residuals (use reference value where masked)
        resid = cumDisp_looked + insarM;          % metres  (may be NaN)
        meanResMat(i,j) = mean(resid,'omitnan');
        residPlot = resid;                        % copy for plotting
        idxMask   = isnan(insarM);                % stations removed by mask
        residPlot(idxMask) = cumDisp_looked(idxMask) + insarRef(idxMask);

        % 5. Plot in the appropriate subplot
        subplot(nrows,ncols,(i-1)*ncols + j);  cla; hold on; box on

        % un-masked first (black), then masked (red)
        plot(r(~idxMask), residPlot(~idxMask), 'k.', 'MarkerSize',18);
        plot(r(idxMask),  residPlot(idxMask),  'r.', 'MarkerSize',18);

        % annotate each point
        for k = 1:nStat
            txtCol = 'k'; if idxMask(k), txtCol = 'r'; end
            text(r(k)+0.02*range(r), residPlot(k), GPSNameList{k}, ...
                 'Color',txtCol, 'FontSize',11, ...
                 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
        end

        xlabel('Distance from caldera center (km)');
        ylabel('GPS – InSAR disp. (m)');
        title(sprintf('coh ≥ %.2f   NaN ≤ %.0f%%', ...
              coh_list(i), nan_list(j)*100), 'FontSize',11);
        grid on
    end
end
sgtitle('Residual plots — masked GPS stations in red','FontSize',16);


%% summary heat-map of mean residual
figure(21); clf;
imagesc(nan_list*100,  coh_list,  flipud(abs(meanResMat)));
set(gca,'YDir','normal');
clim([0.17, 0.226]);
colorbar;
xlabel('Allowed NaN-pixels per block  (%)');
ylabel('Coherence threshold');
title('Mean GPS – InSAR residual (m)');


%% Test coherence thresholds / methods
figure(5);
coh_list = [0.2, 0.3, 0.4];
coh_names = ["Mean coherence", "Min coherence", "Bottom 30% coherence"];
nrows = length(coh_list);
ncols = 3;

clf;
for i = 1:nrows
    for j = 1:ncols
        subplot(nrows, ncols, (i-1)*ncols + j);  % Proper subplot indexing
        filtered_resampled_test = downsample_coh(coh_list(i), squeeze(cc_tests(j,:,:)), resampled, S);
        imagesc(longitudes, latitudes, filtered_resampled_test);
        set(gca, 'YDir', 'normal');
        hold on;
        colormap default
        plot(coast_new(:, 1)', coast_new(:, 2)', 'k.', 'HandleVisibility','off');
        title("Threshold = " + coh_list(i) + ", " + coh_names(j));
        xlim([-155.33, -155.22]);
        ylim([19.34, 19.48]);
    end
end


% axis equal;
% colormap gray;

%% Export data
% Initialize an output array
data_out = [];  % Each row will be: [longitude, latitude, displacement, blockSize]

% Get the unique block sizes from the quadtree structure (process largest blocks first)
blockSizes = full(unique(S(S > 0)));
blockSizes = sort(blockSizes, 'descend');
blocks = [];
insarx = [];
insary = [];

for bs = blockSizes'
    % Find the upper-left corners of blocks of the current size
    [r, c] = find(S == bs);
    for k = 1:length(r)
        % Compute the center index for this block
        % (Use floor(bs/2) to get an index offset; if bs is even, this picks the lower center)
        row_center = r(k) + floor(bs/2);
        col_center = c(k) + floor(bs/2);
        
        % Retrieve the displacement for the block (it should be the same for the whole block)
        disp_val = filtered_resampled(row_center, col_center);
        
        % Optionally, skip blocks that were filtered out (if they are NaN)
        if isnan(disp_val)
            continue;
        end
        
        % Get the corresponding geographic coordinates.
        % If longitudes and latitudes are vectors:
        x = longitudes(col_center);   % columns map to x (longitude)
        y = latitudes(row_center);      % rows map to y (latitude)
        
        % Append the data: [x, y, displacement, blockSize]
        data_out = [data_out; x, y, disp_val, bs];
        insarx = [insarx, x];
        insary = [insary, y];
        blocks = [blocks, bs];
    end
end

% Save the data to a text file (or CSV)
writematrix(data_out, 'ascending.txt', 'Delimiter', 'tab');
fprintf('Saved %d quadtree resampled points with geographic coordinates and look vector displacements.\n', size(data_out,1));
%% Export data for cov matrix generation
cum_disp_rad = cum_disp_asc*2*pi/ lambda;
save asc_cov_data.mat X Y latitudes longitudes avg_cc_asc cum_disp_rad blocks insarx insary

%% Displacement video 
doDispVideo = 0; 
if doDispVideo == 1
v = VideoWriter('kilauea_timeseries.mp4','MPEG-4');
v.FrameRate = 6;
open(v);
figure;
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf, 'Position',[160 440 664 584]);
set(gcf, 'Position',[160 440 864 584]);
for k = 1:nigrams
    cmax = 0.5*max(max(max(abs(displacement))));
    rgb = dishgt_nocycle(disp_amp(:,:,k).',displacement(:,:,k).',cmax,naz,nr);
    imagesc(longitudes, latitudes, rgb);
    title(strcat("Cumulative Displacement, ",string(t_date(k))," to " ,string(t_date(k+1))));
    disp_cbar = colorbar; ylabel(disp_cbar,'Displacement (cm)','FontSize',12,'FontName','Segoe UI','Rotation',90);
    colormap(diverging_colormap); 
    clim([-cmax cmax]);
    set(gca,'YDir','normal')
    ax = gca; ax.FontName = 'Segoe UI'; ax.FontSize = 14;
    axis image;
    xlim([-155.35 -155.1]); ylim([19.22, 19.5]);

    F = getframe(gcf);
    writeVideo(v,F);
    % pause(0.2)
end
close(v);

end

%% Displacement median filter 
medfilt_disp = zeros(size(displacement));
for k = 1:nigrams
    medfilt_disp(:,:,k) = medfilt2(displacement(:,:,k), [50 50]);
end
%% Video of phase and displacement 
v = VideoWriter('kilauea_both_timeseries.mp4','MPEG-4');
v.FrameRate = 2;
open(v);
figure;
% set(gcf, 'Position', get(0, 'Screensize'));
set(gcf, 'Position',[160 440 1064 584]);
for k = 1:nigrams
    subplot(1,2,1)
    rgb = dismph(acc, angle(wrapped_int(:,:,k)), pi,naz, nr);
    imagesc(longitudes, latitudes, rgb);
    title(strcat("Phase Difference, ",string(t_date(k))," to " ,string(t_date(k+1))));
    clim([-pi pi]);
    set(gca,'YDir','normal')
    axis image; 
    ax = gca; ax.FontName = 'Segoe UI'; ax.FontSize = 14;
    phase_cbar = colorbar; ylabel(phase_cbar,'Phase (rad)','FontSize',12,'FontName','Segoe UI','Rotation',90);
    colormap(ax(1));
    xlim([-155.35 -155.1]); ylim([19.22, 19.5]);

    subplot(1,2,2)
    cmax = 0.5*max(max(max(abs(displacement))));
    rgb = dishgt(disp_amp(:,:,k).',displacement(:,:,k).',cmax,naz,nr);
    imagesc(longitudes, latitudes, rgb);
    title(strcat("Cumulative Displacement, ",string(t_date(k))," to " ,string(t_date(k+1))));
    disp_cbar = colorbar; ylabel(disp_cbar,'Displacement (cm)','FontSize',12,'FontName','Segoe UI','Rotation',90);
    clim([-cmax cmax]);
    set(gca,'YDir','normal')
    ax = gca; ax.FontName = 'Segoe UI'; ax.FontSize = 14;
    colormap(ax(1), flipud(diverging_colormap));
    hold on; 
    contour(longitudes, latitudes, medfilt_disp(:,:,k).',[-80 -70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80 100 150 200 250 300],'k','LineWidth',1);
    axis image; 
    xlim([-155.35 -155.1]); ylim([19.22, 19.5]);


    F = getframe(gcf);
    writeVideo(v,F);
    % pause(0.2)
end
close(v);


%% Display in a way you can see values by clicking on the image 
% Change which image it is by changing (:,:,1) or (:,:,end) to other points
% besides the first and last in the time series. 
read_unw = 1; 
if read_unw
    figure; imagesc(longitudes, latitudes, displacement(:,:,end).');
    colorbar;
    title('Unwrapped Phase, May-August')
    set(gca,'YDir','normal')

    figure; imagesc(longitudes,latitudes,cc_array(:,:,1));
    title('Coherence of First Interferogram (early May)')
    colorbar; colormap gray;
    set(gca,'YDir','normal')

    figure; imagesc(longitudes,latitudes, angle(wrapped_int(:,:,1)));
    title('Wrapped phase of First Interferogram (early May)')
    colorbar; colormap turbo
    set(gca,'YDir','normal')
end


%% function for nice display 
function rgb = dishgt(amp, phase, cmax, nr, naz)
    upramp = linspace(100,255,120)/255;
    dnramp = linspace(255,100,120)/255;
    oneramp = ones(120,1);
    zeroramp = zeros(120,1);
    r = zeros(360,1);
    g = zeros(360,1);
    b = zeros(360,1);
    r(1:120)=upramp; 
    r(121:240)=oneramp; 
    r(241:360)=dnramp;
    g(1:120)=dnramp; 
    g(121:240)=upramp; 
    g(241:360)=oneramp;
    b(1:120)=oneramp; 
    b(121:240)=dnramp; 
    b(241:360)=upramp;
    % fig = plt.figure(figsize=(10, 10)) 
    % phase = angle(image);
    phase = phase*180/cmax; % originally *180/pi for dismph
    exponent=0.3;
    scale = 1;
    mag=abs(amp).^exponent; 
    ampsum=sum(sum(mag));
    ampi=naz*nr;
    % ampi=sum(sum(mag~=0));
    scalemag = (scale*150/(ampsum/ampi))/256;
    % mag = np.clip(mag*scalemag,0,1);
    mag = mag*scalemag; 
    mag(mag>1) = 1; 
    mag(mag<0) = 0;
    icolor = round(phase); %phase.astype(int);
    % icolor = np.where(icolor < 0, icolor+360, icolor);
    icolor(icolor>360) = 360;
    icolor(icolor<1) = 1;
    % load color table values for image
    rgb = zeros(nr,naz,3) ;
    rgb(:,:,1) = (r(icolor(:,:)).*mag);
    rgb(:,:,2) = (g(icolor(:,:)).*mag);
    rgb(:,:,3) = (b(icolor(:,:)).*mag);
    % figure; 
    % imagesc(rgb); axis image;
end


%% function for nice display 
function rgb = dismph(amp, phase, cmax, nr, naz)
    upramp = linspace(100,255,120)/255;
    dnramp = linspace(255,100,120)/255;
    oneramp = ones(120,1);
    zeroramp = zeros(120,1);
    r = zeros(360,1);
    g = zeros(360,1);
    b = zeros(360,1);
    r(1:120)=upramp; 
    r(121:240)=oneramp; 
    r(241:360)=dnramp;
    g(1:120)=dnramp; 
    g(121:240)=upramp; 
    g(241:360)=oneramp;
    b(1:120)=oneramp; 
    b(121:240)=dnramp; 
    b(241:360)=upramp;
    % fig = plt.figure(figsize=(10, 10)) 
    % phase = angle(image);
    phase = phase*180/cmax; % originally *180/pi for dismph
    min(min(phase))
    exponent=0.3;
    scale = 1;
    mag=abs(amp).^exponent; 
    ampsum=sum(sum(mag));
    ampi=naz*nr;
    % ampi=sum(sum(mag~=0));
    scalemag = (scale*150/(ampsum/ampi))/256;
    % mag = np.clip(mag*scalemag,0,1);
    mag = mag*scalemag; 
    mag(mag>1) = 1; 
    mag(mag<0) = 0;
    icolor = round(phase + 180); %phase.astype(int);
    % icolor = np.where(icolor < 0, icolor+360, icolor);
    icolor(icolor>360) = 360;
    icolor(icolor<1) = 1;
    % load color table values for image
    rgb = zeros(nr,naz,3) ;
    rgb(:,:,1) = (r(icolor(:,:)).*mag);
    rgb(:,:,2) = (g(icolor(:,:)).*mag);
    rgb(:,:,3) = (b(icolor(:,:)).*mag);
    % figure; 
    % imagesc(rgb); axis image;
end
