%% This script converts the formatting from the output of the MCMC analysis to the spheroid.m formatting and vice versa

function m = get_full_m(prior, opt_params, forward, pressure_type)
% Prior - values to fill in the blanks of optimization result
% opt_params - MCMC output
% forward = true means that we want to get a full 1x16 m state vector to
% input into the spheroid.m code based on our somewhat incomplete
% optimization result
% forward = false means that we want to get the formatting of the output
% MCMC vector from an input 1x16 m state vector.
% pressure_type - "insar" or "gps". Tells which pressure change to use for
% the output m vector.
%
% Variable names: ["dpHMM_insar", "dpHMM_gps", "volHMM", 'xHMM', 'yHMM', 'zHMM', 'HMM aspect ratio"
% "xSC', ySC, dSC, "SC aspect ratio", "dip", strike, "dpSC_insar", "dpSC_gps", "volSC"];

    % Get full m array based on optimized params
    if(forward)
        % Get geometry based on insar pressure change
        if(pressure_type == "insar")
            opt_params = [opt_params(1), opt_params(3:14), opt_params(16)];%, opt_params(15)];
            % opt_params = [opt_params(1), opt_params(3:end-2), opt_params(end)];
        % Otherwise we want the GPS pressure change
        else
            % opt_params = [opt_params(2:end-3), opt_params(end-1:end)];
            opt_params = [opt_params(2), opt_params(3:12), opt_params(15:16)];%, opt_params(15)];
        end
        
        % Fixed aspect ratio for HMM
        aspect_ratio_HMM = opt_params(6); %1.7496;
        % Setting parameters based on combination of prior and MCMC soln
        HMM_volume = opt_params(2);
        opt_vert_sd = (3/(4*pi) * HMM_volume * (aspect_ratio_HMM^2))^(1/3);
        opt_horiz_sd = opt_vert_sd/(aspect_ratio_HMM);
        % Output HMM geometry vector
        mHMM = [opt_vert_sd, opt_horiz_sd, prior(3:4), opt_params(3), opt_params(4), ...
            opt_params(5) - abs(opt_vert_sd), opt_params(1)];
        
        % Same process for SC but fewer parameters are retrieved from the prior
        aspect_ratio_SC = opt_params(10);
        SC_volume = opt_params(end);
        opt_vert_sd = (3/(4*pi) * SC_volume * (aspect_ratio_SC^2))^(1/3);
        opt_horiz_sd = opt_vert_sd/(aspect_ratio_SC);
        prior_SC = prior(10:end);
        mSC = [opt_vert_sd, opt_horiz_sd, opt_params(11:12), opt_params(7:9), opt_params(end-1)];
        
        % Full output geometry vector containing both HMM and SC
        m = real([mHMM, mSC]);
        
    % Get optimization result vector from an inputted m vector
    else
        m = real([prior(8), prior(8), (4/3) * pi * prior(1) * prior(2)^2, prior(5), prior(6), prior(7) + prior(1), prior(1)/prior(2), ...
            prior(13:15), prior(9)/prior(10), prior(11:12), prior(end), prior(end), (4/3) * pi * prior(9) * prior(10)^2]);
    end

end