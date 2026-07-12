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

fid=fopen('20180508_20180806.cc','r');
dat=fread(fid,[nr*2 naz],'float32');
fclose(fid);
cc = dat(nr+1:end, :).'; % InSAR Coherence 

fid=fopen('20180508_20180806.unw','r');
dat=fread(fid,[nr*2 naz],'float32');
fclose(fid);
unw = dat(nr+1:end, :).'; % Unwrapped phase in radians (check units!)

acc=dat(1:nr,:).';
meanacc=mean(acc(:));
stdacc=std(acc(:));
acc(acc>meanacc+1*stdacc)=meanacc+1*stdacc;
acc(acc<meanacc-1*stdacc)=meanacc-1*stdacc;
acc=acc-(max(0,meanacc-3*stdacc));
acc=acc/max(acc(:));
figure;
imagesc(acc'); colormap gray;
axis image;

fprintf('Interferometric Data loaded.\n');
%% One option to display in a nice-looking way
rgb = dishgt(acc, unw, max(max(unw)),naz, nr);
figure; imagesc(longitudes, latitudes, rgb); 
title('Unwrapped interferogram')
set(gca,'YDir','normal')

rgb = dismph(acc, angle(wrapped_int), pi,naz, nr);
figure; imagesc(longitudes, latitudes, rgb); 
title('Wrapped interferogram')
set(gca,'YDir','normal')

%% Display in a way you can see values by clicking on the image 
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