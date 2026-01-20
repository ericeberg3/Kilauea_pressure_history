%% Creates tilt greens functions

function [gTiltHMM, gTiltSC] = createtiltgreens(mHMM, mSC, dtheta, keep, mu, flag)
    xyll = [-155.2937; 19.3905];
    xy = llh2local(xyll, [-155.2784, 19.4073]) * 1000;
    x = xy(1);
    y = xy(2);
    z = zeros(size(x));

    if nargin < 5
        mu = 3.08*10^9;
    end
    if nargin < 6
        flag = 'pressure';
    end
    % Sign convention negative dp is positive pressure increase
    
    npitloc = coord('NPIT', 'llh');
    npitloc = llh2local(npitloc(1:2), [-155.2784, 19.4073]) * 1000;
    
    [~, dHMM, ~, ~] = spheroid(mHMM, [x; y; z], 0.25, 3.08*10^9, flag);
    [~, dSC, ~, ~] = spheroid(mSC, [x; y; z], 0.25, 3.08*10^9, flag);
    dHMM = real(dHMM);
    dSC = real(dSC);
    
    % Convert rad to urad
    % Negative sign included to put the atan in the correct quadrant
    gTiltHMM = -real([atan2(dHMM(3), 1), atan2(dHMM(6), 1)]) .* 1e6; % real([atan2(dHMM(3) * cos(dtheta) - dHMM(6) * sin(dtheta), 1), atan2(dHMM(3) * sin(dtheta) + dHMM(6) * cos(dtheta), 1)]) .* 1e6;
    gTiltSC = -real([atan2(dSC(3), 1), atan2(dSC(6), 1)]) .* 1e6; % real([atan2(dSC(3) * cos(dtheta) - dSC(6) * sin(dtheta), 1), atan2(dSC(3) * sin(dtheta) + dSC(6) * cos(dtheta), 1)]) .* 1e6;
    
    
    if(keep)
        save Data/gTiltHMM.mat
        save Data/gTiltSC.mat
    end
end