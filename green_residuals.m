function res = green_residuals(m, m_guess, x, y, z, u, recstds, tilt, nanstatsbeg, priormeans, priorvariances)

npitloc = coord('NPIT', 'llh');
npitloc = llh2local(npitloc(1:2), [-155.2784, 19.4073]) * 1000;
mHMM = [m_guess(1:7), m(1)];
% mHMM = [1600.79, 914.47, 90, 0, 70.3, 183, -1940, m(1)];
mSCguess = m_guess(9:end); % [277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, 1e7];
mSC = [m(2:4), mSCguess(4:end-1), m(end)];
% consts = m(17:end);

% tilt(1) = tilt(1) * cos(tiltoffset(1)) - sin(tiltoffset(1));
% tilt(2) = tilt(2) * cos(tiltoffset(1)) + sin(tiltoffset(1));

constsfull = zeros(length(nanstatsbeg), 3);
j = 1;
for i = 1:length(nanstatsbeg)
    if(nanstatsbeg(i) == 1)
        % constsfull(i, :) = consts(j);
        j = j + 1;
    end
end
clear j i consts

% Forward model computations
[gHMM, ~, ~, ~] = spheroid(mHMM, [x(1:end-1); y(1:end-1); z(1:end-1)], 0.25, 3.08*10^9);
[gSC, ~, ~, ~] = spheroid(mSC, [x(1:end-1); y(1:end-1); z(1:end-1)], 0.25, 3.08*10^9);
[~, dHMM, ~, ~] = spheroid(mHMM, [x(end); y(end); z(end)], 0.25, 3.08*10^9);
[~, dSC, ~, ~] = spheroid(mSC, [x(end); y(end); z(end)], 0.25, 3.08*10^9);

% Return large residual if nan value
if(~isreal(dHMM) || ~isreal(dSC))
    res = 1e10;
    return;
end

gtot = gHMM + gSC;
[gTiltHMM, gTiltSC] = createtiltgreens(mHMM, mSC, 0, false); % [atan2(real(dHMM(3)), 1), atan2(real(dHMM(6)), 1)] .* 1e6 + [atan2(real(dSC(3)), 1), atan2(real(dSC(6)), 1)] .* 1e6;
gTilt = gTiltHMM + gTiltSC;
rest = tilt - gTilt;
resd = u + constsfull - gtot';

% Compute total residual 
weights = [recstds(1) .* ones(1, length(x) - 1), recstds(2) .* ones(1, length(x) -1), recstds(3) .* ones(1, length(x) - 1), recstds(4:5)];
% weights = [recstds(1) .* ones(1, length(x) - 1), recstds(2) .* ones(1, length(x) -1), recstds(3) .* ones(1, length(x) - 1)];
% weights(end-1:end) = 0; Include line if we want to not include tilt
W = diag(weights);
res = ([reshape(resd.',1,[]), 0 * reshape(rest.',1,[])]); % NOTE: tilt residual is NOT included due to long term drift concerns

% Gaussian prior penalty
prior_residuals = (m - priormeans).^2 ./ priormeans.^2;
prior_penalty = sum(prior_residuals);

res = real(res * (W.^2 * res') * 1e-5); %  + prior_penalty;
% res = prior_penalty;

% disp("Sim tilt: " + gTilt + " Real tilt: " + tilt + "\n")
% disp(dot(res, res'))

end