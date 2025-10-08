%% This function calculates the pressure at all timesteps using a single large least squares step

function dp = TimeDependentLSQtilt(gHMMflat, gSCflat, gTiltHMM, gTiltSC, ux, uy, uz, tiltx, tilty, recstds, ...
    GPSNameList, rwsigma, dp_weight, useBoundaryRows)
    gtotHMM = cat(1, gHMMflat);
    gtotSC = cat(1, gSCflat);
    % Getting rid of nan stations
    nanstations = [isnan(ux(:, 1200)); isnan(ux(:, 1200)); isnan(ux(:, 1200)); isnan(tiltx(1200)); isnan(tilty(1200))];
    nanstatbeginning = [isnan(ux(:, 1)); isnan(uy(:, 1)); isnan(uz(:, 1))];
    nanstations(find(contains(GPSNameList,'UWEV')):length(GPSNameList):end) = 1;
    nanstations(find(contains(GPSNameList,'BYRL')):length(GPSNameList):end) = 1;
    nanstations(find(contains(GPSNameList,'CRIM')):length(GPSNameList):end) = 1;
    nanstations = logical(nanstations);

    gtotHMM(nanstations) = [];
    gtotSC(nanstations) = [];
    nstat = size(ux, 1);
    ux(nanstations(1:nstat), :) = [];
    uy(nanstations(nstat+1:nstat*2), :) = [];
    uz(nanstations(nstat*2+1:nstat*3), :) = [];
    GPSNameList = GPSNameList(~nanstations(1:length(GPSNameList)));

    %% Set up weight matrix  & green's functions matrix G
    % Wdiag = the diagonal of the GPS components of the Sigma^-1 matrix
    % tiltW = the tilt block of the Sigma^-1 matrix
    noffset = sum(nanstatbeginning);
    ntime = length(ux);
    greensize = length(gtotHMM);
    tiltsize = length(gTiltHMM);

    % Expand G if imposing the offset BC
    if(useBoundaryRows)                                                
        G = zeros((greensize + 2) * ntime + 2, ntime * 2 + noffset);
        ureal = zeros(1, (greensize + 2) * ntime + 2);
    else
        G = zeros((greensize + 2) * ntime, ntime * 2 + noffset); 
        ureal = zeros(1, (greensize + 2) * ntime);
    end
    
    Wdiag = ones(size(greensize*ntime, 2), 1);
    tiltW = zeros(tiltsize*ntime, tiltsize*ntime);

    %% Set up offsets column of G
    offset_matrix = zeros(length(G), noffset);
    naninds = find(nanstatbeginning);
    for i = 1:length(find(nanstatbeginning))
        offset_matrix(naninds(i):length(nanstatbeginning):end, i) = 1;
    end
    % Fill in G matrix with offsets section (according to the order of
    % GPSNameList)
    G(:, (ntime*2 + 1):end) = offset_matrix;

    %% Fill in G and W matrices 

    for i = 1:ntime
        % Fill in GPS stations greens functions into G and weights into
        % Wdiag
        G((((i-1) * greensize) + 1):(i * greensize), i) = gtotHMM';
        G(((i-1) * greensize + 1):(i * greensize), i + ntime) = gtotSC';
        ureal((((i -1) * greensize) + 1):((i - 1) * greensize + greensize)) = cat(1, ux(:, i), uy(:, i), uz(:, i));
        sizegps = greensize;
        Wdiag((greensize*(i - 1) + 1) : (greensize*(i - 1) + sizegps/3)) = recstds(end - 3*(length(recstds) - 2)/3); % weight x component
        Wdiag((greensize*(i - 1) + 1 + sizegps/3) : (greensize*(i - 1) + 2*sizegps/3)) = recstds(end - 2*(length(recstds) - 2)/3); % weight y component
        Wdiag((greensize*(i - 1) + 1 + 2*sizegps/3) : (greensize*(i - 1) + 3*sizegps/3)) = recstds(end - 1*(length(recstds) - 2)/3); % weight z component

        % Fill in tilt stations greens functions into G and weights into
        % tiltW
        gpsendind = greensize*ntime + 1;
        % G((gpsendind + tiltsize*(i-1)):(gpsendind + tiltsize*(i-1) + 1), i) = gTiltHMM';
        % G((gpsendind + tiltsize*(i-1)):(gpsendind + tiltsize*(i-1) + 1), i + ntime) = gTiltSC';
        G((gpsendind + i-1), i) = gTiltHMM(1);
        G((gpsendind + length(tiltx) + i - 1), i) = gTiltHMM(2);
        G((gpsendind + i-1), i + ntime) = gTiltSC(1);
        G((gpsendind + length(tiltx) + i-1), i + ntime) = gTiltSC(2);
        
        ureal(gpsendind + i-1) = tiltx(i);
        ureal(gpsendind + length(tiltx) + i-1) = tilty(i);
        % ureal((gpsendind + tiltsize*(i-1)):(gpsendind + tiltsize*(i-1) + 1)) = [tiltx(i), tilty(i)];
        deltat = 1; % in hours
        % The random walk standard deviation sigma that I am using is the
        % standard deviation from the no-signal period (urad/sqrt(hr))
        % Sigma_{i, i}
        tiltW(i, i) = 1/(2*rwsigma^2) * (2/deltat);
        tiltW(i + ntime, i+ntime) = 1/(2*rwsigma^2) * (2/deltat);

        % % Sigma_{i - 1, i}
        tiltW(i -1 + ntime, i + ntime) = -1/(2*rwsigma^2) * (1/deltat);
        if(i > 1)
            tiltW(i-1, i) = -1/(2*rwsigma^2) * (1/deltat);
        end

        % Sigma_{i + 1, i}
        tiltW(i + 1, i) = -1/(2*rwsigma^2) * (1/deltat);
        if(i < ntime)
            tiltW(i + ntime + 1, i + ntime) = -1/(2*rwsigma^2) * (1/deltat);
        end
    end

    % Fill in the two 0 pressure conditions into G
    if(useBoundaryRows)
        G(end-1, 1) = 1e3;
        G(end, ntime + 1) = 1e3;
    end

    % Set specific boundary cases of RW weighting
    % disp("RW Sigma = " + rwsigma)
    % make the random walk "start over" for the next component of tilt
    tiltW(ntime, ntime + 1) = 0;
    tiltW(ntime + 1, ntime) = 0;
    tiltW(1, 1) = 1/(2*rwsigma^2) * (1/deltat);
    tiltW(ntime, ntime) = 1/(2*rwsigma^2) * (1/deltat);
    tiltW(ntime + 1, ntime+1) = 1/(2*rwsigma^2) * (1/deltat);
    tiltW(2*ntime, 2*ntime) = 1/(2*rwsigma^2) * (1/deltat);
    
    n = size(ureal, 2);
    len_Wdiag = length(Wdiag);
    
    % Construct the diagonal block with Wdiag
    diag_indices = (1:len_Wdiag)';
    Wdiag_values = Wdiag(:);
    [tilt_rows, tilt_cols, tilt_vals] = find(tiltW);
    
    % Adjust indices to position tiltW in W
    tilt_rows = tilt_rows + len_Wdiag;
    tilt_cols = tilt_cols + len_Wdiag;
    
    % Combine indices and values and include dp condition weighting
    
    W_rows = [diag_indices; tilt_rows];
    W_cols = [diag_indices; tilt_cols];
    W_vals = [Wdiag_values; tilt_vals];
    
    % Create W as a sparse matrix
    
    W = sparse(W_rows, W_cols, W_vals, n, n);
    if(useBoundaryRows)
        W(end-1, end-1) = dp_weight;
        W(end, end) = dp_weight;
    end
    
    % Convert G to sparse representation
    G = sparse(G);
    
    % Transpose ureal if needed
    ureal = ureal';

    % Clean up NaN values in preparation of least squares inversion
    G(isnan(ureal), :) = [];
    W(isnan(ureal), :) = [];
    W(:, isnan(ureal)) = [];
    ureal(isnan(ureal)) = [];

    % Doing the weighted least squares
    A = (G' * W * G);
    dp = (A \ (G' * W * ureal))';
end