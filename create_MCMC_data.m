function data = create_MCMC_data(m, m_guess, x, y, z, insarx, insary, ...
    look, insar_lengths, nanstatsbeg)

    npitloc = coord('NPIT', 'llh');
    npitloc = llh2local(npitloc(1:2), [-155.2784, 19.4073]) * 1000;
    
    % aspect_ratio = 1.7496;
    % opt_vert_sd = (3/(4*pi) * m(1) * (aspect_ratio^2))^(1/3);
    % opt_horiz_sd = opt_vert_sd/(aspect_ratio);
    
    % Create m based on the gps pressure change
    m_full = get_full_m(m_guess, m, true, "gps");
    mHMM = m_full(1:8);
    mSC = m_full(9:end);
    
    % tilt(1) = tilt(1) * cos(tiltoffset(1)) - sin(tiltoffset(1));
    % tilt(2) = tilt(2) * cos(tiltoffset(1)) + sin(tiltoffset(1));
    
    j = 1;
    for i = 1:length(nanstatsbeg)
        if(nanstatsbeg(i) == 1)
            % constsfull(i, :) = consts(j);
            j = j + 1;
        end
    end
    clear j i consts
    
    % Forward model computations
    % GPS data generation
    [gHMM, ~, ~, ~] = spheroid(mHMM, [x(1:end); y(1:end); z(1:end)], 0.25, 3.08*10^9, 'volume');
    [gSC, ~, ~, ~] = spheroid(mSC, [x(1:end); y(1:end); z(1:end)], 0.25, 3.08*10^9, 'volume');
    
    gtot = gHMM + gSC;
    GPS_data = gtot(:);

    % Create m based on the insar pressure change
    m_full = get_full_m(m_guess, m, true, "insar");
    mHMM = m_full(1:8);
    mSC = m_full(9:end);

    insar_data = zeros(size(insarx));
    
    % InSAR data generation
    for i = 1:length(insar_lengths)
        if i == 1; inds = 1:insar_lengths(i); end
        if i == 2; inds = (insar_lengths(i-1) + 1):sum(insar_lengths); end

        [gHMM, ~, ~, ~] = spheroid(mHMM, [insarx(inds); insary(inds); zeros(size(insarx(inds)))], 0.25, 3.08*10^9, 'volume');
        [gSC, ~, ~, ~] = spheroid(mSC, [insarx(inds); insary(inds); zeros(size(insarx(inds)))], 0.25, 3.08*10^9, 'volume');
    
        gtot = gHMM + gSC;
        insar_data(inds) = gtot' * look(:,i);
    end
    
    data = real([GPS_data; insar_data(:)]);
end