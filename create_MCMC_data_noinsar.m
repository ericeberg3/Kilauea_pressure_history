function data = create_MCMC_data_noinsar(m, m_guess, x, y, z, recstds, nanstatsbeg)

npitloc = coord('NPIT', 'llh');
npitloc = llh2local(npitloc(1:2), [-155.2784, 19.4073]) * 1000;

aspect_ratio = 1.7496;
opt_vert_sd = (3/(4*pi) * m(1) * (aspect_ratio^2))^(1/3);
opt_horiz_sd = opt_vert_sd/(aspect_ratio);

mHMM = [opt_vert_sd, opt_horiz_sd, m_guess(3:6), m_guess(7) - abs(opt_vert_sd - m_guess(1)), m(2)];
% mHMM = [1600.79, 914.47, 90, 0, 70.3, 183, -1940, m(1)];
mSCguess = m_guess(9:end); % [277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, 1e7];
mSC = [m(2:4), mSCguess(4:end-1), m(end)];
% consts = m(17:end);

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
[gHMM, ~, ~, ~] = spheroid(mHMM, [x; y; z], 0.25, 3.08*10^9);
[gSC, ~, ~, ~] = spheroid(mSC, [x; y; z], 0.25, 3.08*10^9);

gtot = gHMM + gSC;
data = gtot';
end