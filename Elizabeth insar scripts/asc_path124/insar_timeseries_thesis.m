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
%%
% cc_array = zeros(naz,nr,nigrams);
% for k = 1:nigrams
%     fprintf(strcat('Loading phases_', num2str(k),'...\n'));
%     fid=fopen(strcat(intlist{k}(1:17),'.cc'),'r');
%     phase_cc=fread(fid,[inf],'float32');
%     phase_cc = reshape(phase_cc, [2*nr naz]); 
%     % clear phase_cc; % clears up space
%     cc_array(:,:,k)=phase_cc(nr+1:end,:).'; % Wrapped interferogram
%     fclose(fid);
% end
%%
read_cc_unw = 0;
if read_cc_unw
    fid=fopen('20180508_20180806.cc','r');
    dat=fread(fid,[nr*2 naz],'float32');
    fclose(fid);
    cc = dat(nr+1:end, :).'; % InSAR Coherence
end

if read_cc_unw
    fid=fopen('20180508_20180806.unw','r');
    dat=fread(fid,[nr*2 naz],'float32');
    fclose(fid);
    unw = dat(nr+1:end, :).'; % Unwrapped phase in radians (check units!)
end

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
%% One option to display in a nice-looking way
if read_cc_unw
    rgb = dishgt(acc, unw, max(max(unw)),naz, nr);
    figure; imagesc(longitudes, latitudes, rgb);
    title('Unwrapped interferogram')
    set(gca,'YDir','normal')
end

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
doPhaseVideo = 0; 
if doPhaseVideo == 1
v = VideoWriter('kilauea_timeseries.mp4','MPEG-4');
v.FrameRate = 6;
open(v);
figure;
% set(gcf, 'Position', get(0, 'Screensize'));
set(gcf, 'Position',[160 440 664 584]);
for k = 1:nigrams
    rgb = dismph(acc, angle(wrapped_int(:,:,k)), pi,naz, nr);
    imagesc(longitudes, latitudes, rgb);
    title(strcat("Phase Difference, ",string(t_date(k))," to " ,string(t_date(k+1))));
    colorbar; colormap(howard_colormap); 
    clim([-pi pi]);
    set(gca,'YDir','normal')
    axis image; 
    xlim([-155.35 -155.1]); ylim([19.22, 19.5]);

    F = getframe(gcf);
    writeVideo(v,F);
    % pause(0.2)
end
close(v);
end

%% Load displacement
fprintf(strcat('Loading Displacement...\n'));
fid=fopen("displacement",'r');
dat=fread(fid,[inf],'float32');
fclose(fid);
disp=reshape(dat,nr*2,naz,nigrams); % Displacement
disp_amp = disp(1:nr,:,:);
displacement = disp(nr+1:end,:,:);

lambda = 5.56;
displacement = displacement/(4*pi)*lambda; 

% Mask 
[X,Y] = meshgrid(1:size(disp_amp(:,:,end),1),1:size(disp_amp(:,:,end),2));
mask = 1020*X - 209*Y < 586291; 
% figure; imagesc(disp_amp(:,:,end) + disp_amp(:,:,end).*mask.');
% clim([0 1])
disp_amp = disp_amp + disp_amp .* mask.';

cmax = 0.5*max(max(max(abs(displacement))));
rgb = dishgt_nocycle(disp_amp(:,:,end).',displacement(:,:,end).',cmax,naz,nr);
figure; imagesc(longitudes, latitudes, rgb);
set(gca,'YDir','normal')
colormap(diverging_colormap); colorbar; clim([-cmax cmax]); 
% figure; imagesc(disp_amp(:,:,1)); 
% figure; imagesc(displacement(:,:,1));
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
    colormap(ax(1), howard_colormap);
    xlim([-155.35 -155.1]); ylim([19.22, 19.5]);

    subplot(1,2,2)
    cmax = 0.5*max(max(max(abs(displacement))));
    rgb = dishgt_nocycle(disp_amp(:,:,k).',displacement(:,:,k).',cmax,naz,nr);
    imagesc(longitudes, latitudes, rgb);
    title(strcat("Cumulative Displacement, ",string(t_date(1))," to " ,string(t_date(k+1))));
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
if read_cc_unw
    figure; imagesc(longitudes, latitudes, unw);
    colorbar;
    title('Unwrapped Phase')
    set(gca,'YDir','normal')

    figure; imagesc(longitudes,latitudes,cc);
    title('Coherence')
    colorbar; colormap gray;
    set(gca,'YDir','normal')

    figure; imagesc(longitudes,latitudes, angle(wrapped_int));
    title('Wrapped phase')
    colorbar; colormap turbo
    set(gca,'YDir','normal')
end

% READ IN DISPLACEMENT PRODUCT! 

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