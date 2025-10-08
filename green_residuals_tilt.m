function res = green_residuals_tilt(dtheta, m, x, y, z, u, recstds, tilt, nanstatsbeg)

mHMM = m(1:8);
mSC = m(9:16);
consts = m(17:end);

% tilt(1) = 1e6 * atan(tan(1e-6 * tilt(1)) * cos(dtheta) - tan(1e-6 * tilt(2)) * sin(dtheta));
% tilt(2) = 1e6 * atan(tan(1e-6 * tilt(2)) * sin(dtheta) + tan(1e-6 * tilt(2)) * cos(dtheta));

constsfull = zeros(length(nanstatsbeg), 3);
j = 1;
for i = 1:length(nanstatsbeg)
    if(nanstatsbeg(i) == 1)
        constsfull(i, :) = consts(j);
        j = j + 1;
    end
end
clear j i consts

[gHMM, ~, ~, ~] = spheroid(mHMM, [x(1:end-1); y(1:end-1); z(1:end-1)], 0.25, 3.08*10^9);
[gSC, ~, ~, ~] = spheroid(mSC, [x(1:end-1); y(1:end-1); z(1:end-1)], 0.25, 3.08*10^9);

[~, dHMM, ~, ~] = spheroid(mHMM, [x(end); y(end); z(end)], 0.25, 3.08*10^9);
[~, dSC, ~, ~] = spheroid(mSC, [x(end); y(end); z(end)], 0.25, 3.08*10^9);
gtot = gHMM + gSC;


gTilt = createtiltgreens(mHMM, mSC, dtheta, false);
gTilt = gTilt(1) + gTilt(2);
rest = tilt - gTilt;

weights = [recstds(1) .* ones(1, length(x) - 1), recstds(2) .* ones(1, length(x) -1), recstds(3) .* ones(1, length(x) - 1), recstds(4:5)];
% weights = [recstds(1) .* ones(1, length(x) - 1), recstds(2) .* ones(1, length(x) -1), recstds(3) .* ones(1, length(x) - 1)];
% weights(end-1:end) = 0; Include line if we want to not include tilt
W = diag(weights);

res = ([reshape(rest.',1,[])]); %used to be 1e-3
% res = ([reshape(resd.',1,[])]);

res = real(res * res') * 1e-4;


% disp("Sim tilt: " + gTilt + " Real tilt: " + tilt + "\n")
% disp(dot(res, res'))

end