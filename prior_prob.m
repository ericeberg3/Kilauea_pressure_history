%%  extractPosteriors_debugFit.m   –– robust baseline & safe ksdensity
clear; clc;

npitloc = coord('NPIT', 'llh');
npitloc = llh2local(npitloc(1:2), [-155.2784, 19.4073]) * 1000;

params = [ ...
    struct('file','Data/post_im/x_hmm_kyle.png',    'name','xHMM',     'xlim', [-100, 200]); ... %[0.30e3 0.50e3] + npitloc(1)); ...
    struct('file','Data/post_im/y_hmm_kyle.png',    'name','yHMM',     'xlim', [0, 400]); ... %[0.20e3 0.50e3] + npitloc(2)); ...
    struct('file','Data/post_im/d_hmm_kyle.png',    'name','dHMM',     'xlim',[-1.3e3, -200]); ... % Test purposes make this from -2e3 to -1e2 [2.05e3 2.20e3] - 1.6e3
    struct('file','Data/post_im/vol_hmm.png',  'name','volHMM',   'xlim',[0 12e9]); ...
    struct('file','Data/post_im/alpha_hmm.png', 'name','alphaHMM','xlim',[1.6 1.8]); ...
    struct('file','Data/post_im/x_sc.png',     'name','xSC',      'xlim',[1.60e3 2.10e3] + npitloc(1)); ...
    struct('file','Data/post_im/y_sc.png',     'name','ySC',      'xlim',[-3.25e3 -2.90e3] + npitloc(2)); ...
    struct('file','Data/post_im/d_sc.png',     'name','dSC',      'xlim',[-3.30e3 -4.20e3]); ...
    struct('file','Data/post_im/strike_sc.png','name','dipSC',    'xlim',[55 70]); ...
    struct('file','Data/post_im/dip_sc.png',   'name','strikeSC', 'xlim',[120 160]); ...
    struct('file','Data/post_im/alpha_sc.png', 'name','alphaSC',  'xlim',[0.10 0.27]); ...
    struct('file','Data/post_im/alpha_sc.png', 'name','dpHMM_insar',  'xlim',[-1e8 1e8]); ...
    struct('file','Data/post_im/alpha_sc.png', 'name','dpHMM_gps',  'xlim',[-1e8 1e8]); ...
    struct('file','Data/post_im/alpha_sc.png', 'name','dpSC_insar',  'xlim',[-1e8 1e8]); ...
    struct('file','Data/post_im/alpha_sc.png', 'name','dpSC_gps',  'xlim',[-1e8 1e8]); ...
    struct('file','Data/post_im/alpha_sc.png', 'name','volSC',  'xlim',[0.5e9 5.5e9]); ...
    struct('file','Data/post_im/mu.png', 'name','mu', 'xlim',[8.5, 11]) % note units of log10(Pa)
];

paramDists = struct();
nP=numel(params); nCol=3; nRow=ceil(nP/nCol);
f=figure('Name','Mask+Fit debug','Color','w','Units','normalized',...
         'Position',[.05 .05 .9 .85]);
tl=tiledlayout(f,nRow,nCol,'TileSpacing','compact');

for p=1:nP
    I=imread(params(p).file); nCols=size(I,2);
    mask=getMask(I); [rows,cols]=find(mask);

    nexttile(tl); imshow(I); hold on;
    if nnz(mask)<10
        title([params(p).name '  – <no pixels>'],'Interpreter','none');
        warning('No pixels for %s',params(p).name); continue
    end
    plot(cols,rows,'.r','MarkerSize',5);

    % ---- height-weighted samples ---------------------------------------
    topRow = accumarray(cols,rows,[nCols 1],@min,NaN);

    darkLine = all(I(:,:,1:3)<40,3);
    darkFrac = sum(darkLine,2)/nCols;
    idxDark  = find(darkFrac>0.05,1,'last');
    idxCurve = max(rows);
    baseline = double( max(idxCurve+2 , ternary(~isempty(idxDark),idxDark,size(I,1))) );

    height = baseline - topRow;  height(height<0|isnan(height))=0;

    edges=linspace(params(p).xlim(1),params(p).xlim(2),nCols+1);
    centres=(edges(1:end-1)+edges(2:end))/2;
    samples=repelem(centres,height');

    if numel(samples)<10
        title(sprintf('%s  – < %d samples>',params(p).name,numel(samples)),...
              'Interpreter','none');
        warning('%s skipped: only %d samples',params(p).name,numel(samples));
        continue
    end

    % ---- histogram ------------------------------------------------------
    if all(samples >= 0)          % every sample is non-negative
        pdfSupp = 'positive';     % 1-sided density
    else
        pdfSupp = 'unbounded';    % allow values on ℝ
    end

    [fhat,xh] = ksdensity(samples, 'Function','pdf', 'Support', pdfSupp);

    plot(x2px(xh,params(p).xlim,nCols), y2px(fhat,baseline,max(height)), ...
         '--','Color',[0 .8 .8],'LineWidth',1.2);

    % ---- fit distribution ----------------------------------------------
    dist=bestFit(samples);
    
    % if strcmp(params(p).name,'alphaHMM')
    %     dist = fitdist(samples','Lognormal');
    %     fprintf('alphaHMM → forced Lognormal\n');
    % end

    
    % >>> special case for volHMM  –– peak-weighted Gamma fit
    if strcmp(params(p).name,'volHMM')
        xdata = samples(:);
    
        % -- step 1: smooth empirical pdf on a fixed grid -------------------
        xGrid = linspace(min(xdata), max(xdata), 120)';
        [fKDE,~] = ksdensity(xdata, xGrid, 'Function','pdf');
    
        % -- step 2: weighted LS objective  (weights = KDE heights) ---------
        w = fKDE.^3 ./ max(fKDE);                      % normalise 0–1
        objFun = @(p) sum( w .* ( gampdf(xGrid,p(1),p(2)) - fKDE ).^2 );
    
        % initial guess from moments (unbiased, but good enough)
        m  = mean(xdata);   v = var(xdata);          % sample mean, var
        k0 = m^2 / v;       th0 = v / m;             % method-of-moments
        p0 = [max(k0,1.2) , th0];                    % ensure k>1 for a mode
    
        % bounds and optimisation
        lb = [1.01 , 1e-6];   ub = [50 , max(xdata)];
        opts = optimset('fminsearch'); opts.MaxIter=4e3; opts.TolX=1e-10;
        pHat = fminsearchbnd(objFun, p0, lb, ub, opts);
    
        k = pHat(1);  th = pHat(2);
    
        dist = makedist('Gamma','a',k,'b',th);      % prob.GammaDistribution
    
        fprintf('volHMM → Gamma  k=%.3f  θ=%.3f  (mode=%.3f)\n', ...
                k,th,(k-1)*th );
    end


    if strcmp(params(p).name,'dpSC_insar') | strcmp(params(p).name,'dpSC_gps')
        a = -1e8;
        b = 1e8;
        dist = makedist('Uniform','lower',a,'upper',b);
    end
    
    % Estimated dpHMM prior to get ~3 MPa in the correct time period
    if strcmp(params(p).name,'dpHMM_insar')
         dist = makedist('Normal','mu',-22e6,'sigma',5e6);
    end
    if strcmp(params(p).name,'dpHMM_gps')
         dist = makedist('Normal','mu',-10e6,'sigma',5e6);
    end

    % Estimate the SC volume prior as a gaussian around 2.5 km^3
    if strcmp(params(p).name,'volSC')
        dist = makedist('Uniform','lower',3.0e9, 'upper', 20e9);
    end
    
    % alphaHMM harder to fit bc of truncation, so manually set mu
    if strcmp(params(p).name,'alphaHMM')
        muHat = 1.78;   % weighted mean
        logL  = @(s) -sum( height' .* log( normpdf(centres, muHat, s) + eps ) );
        sigmaHat = fminbnd(logL, 1e-4, 0.25);
        base = makedist('Normal','mu',muHat,'sigma',sigmaHat);
        dist = truncate(base, params(p).xlim(1), params(p).xlim(2));
    end
    
    xs = linspace(params(p).xlim(1), params(p).xlim(2), 400);
    pf = dist.pdf(xs); scale=max(fhat)/max(pf);
    
    xlabelTxt = '';

    plot(x2px(xs,params(p).xlim,nCols), y2px(pf*scale,baseline,max(height)), ...
         'b','LineWidth',1);
    title(params(p).name,'Interpreter','none'); axis off

    paramDists.(params(p).name)=struct('family',class(dist),'dist',dist,...
        'samples',samples,'xlim',params(p).xlim);
end
save Data/paramDists.mat paramDists
fprintf('\nSaved → paramDists.mat\n');

%% ------------ helpers ----------------------------------------------------
function m = getMask(I)
R = double(I(:,:,1));  G = double(I(:,:,2));  B = double(I(:,:,3));
m = (G - max(R,B)) > 10 & G > 30;        % lower thresholds: pick faint green
end
%% ---------- helper: model selection order (replace old bestFit) ---------
function d = bestFit(s)
cand = {'Lognormal','Gamma','Normal','tLocationScale'};   % skew first!
best = inf;
for c = 1:numel(cand)
    try
        tmp  = fitdist(s', cand{c});
        AIC  = 2*numel(tmp.ParameterValues) - 2*sum(log(pdf(tmp,s')));
        if AIC < best, best = AIC; d = tmp; end
    catch, end
end
if ~exist('d','var')
    d = fitdist(s','Kernel','BandWidth',0.03*range(s));
end
end

function px=x2px(x,xlim,nCols), px=(x-xlim(1))./diff(xlim)*(nCols-1)+.5; end
function py=y2px(y,baseline,maxh), py=baseline - y/max(y)*maxh; end
function o=ternary(c,a,b); if c, o=a; else, o=b; end; end

function y = stableLognormMixPDF(x, gm)
% Vectorised, log-domain evaluation of a log-normal mixture PDF
    x  = x(:);
    lx = log(x);

    % log-pdf of each component:  ln w_k  +  ln φ
    logW  = log(gm.ComponentProportion(:)');
    mu    = gm.mu(:)';          % 1×K
    sig   = sqrt(gm.Sigma(:))'; % 1×K

    logPhi = -0.5*((lx - mu)./sig).^2 - log(sig) - 0.5*log(2*pi);
    logComp = logW + logPhi - lx;          % subtract ln(x) for Jacobian

    % use log-sum-exp for stability
    m  = max(logComp,[],2);
    y  = exp( m + log(sum( exp(logComp - m), 2 )) );
end
% (add once, if not present)
function x = fminsearchbnd(fun,x0,LB,UB,opts)
    x0 = max(LB,min(UB,x0));
    phi=@(y)sin(y); inv=@(y)asin(y);
    x2y=@(x)inv( 2*(x-LB)./(UB-LB) - 1 );
    y2x=@(y)LB + (UB-LB).*(phi(y)+1)/2;
    y0 = x2y(x0);
    y  = fminsearch(@(y)fun(y2x(y)), y0, opts);
    x  = y2x(y);
end


%%
load Data/paramDists.mat
figure(12);
subplot(1,3,1)
x = linspace(1.5, 1.9,600);
plot(x, pdf(paramDists.alphaHMM.dist,x)), title alphaHMM

subplot(1,3,2)
x = linspace(1.5e9,3.5e9,600);
plot(x, pdf(paramDists.volSC.dist,x)),   title volSC

subplot(1,3,3)
x = linspace(-1e8,1e8,600);
stem(x, pdf(paramDists.dpHMM_insar.dist,x)),    title dpHMM (Uniform)

