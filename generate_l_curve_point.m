function [optParams_i, gps_l2, insar_l2, prior_l2] = generate_l_curve_point(i, n_l_curve, u1d, insaru_full, ...
    invStdPWRL, start_params, prior_params, xopt, yopt, zopt, insarx, insary, look, insar_lengths, ...
    nanstatbeginning, cinv_full, gps_weight, insar_weight, prior_weight, paramNames, lb, ub)
    
    disp(['Running Normalized Pattern Search for L-curve point ', num2str(i), '/', num2str(n_l_curve)]);
    
    % 1. Format observed data to match create_MCMC_data output
    data_obs = [u1d(:); insaru_full(:)];
    gps_obs = data_obs(1:length(u1d(:)));
    insar_obs = data_obs(length(u1d(:))+1:end);
    
    % 2. Setup GPS Standard Deviations array
    GPS_std = 1./invStdPWRL;
    sigma = ones(size([u1d(:)]));
    for k = 1:length(u1d(:))
        sigma(k) = GPS_std(mod(k - 1, 3) + 1);
    end
    
    % 3. NORMALIZATION WRAPPER
    % This maps values from [0, 1] back to their true physical bounds [lb, ub]
    norm_to_true = @(m_norm) lb + m_norm .* (ub - lb);
    
    % The objective function now takes the normalized [0,1] array, scales it 
    % to physical bounds, and passes it to your physics solver.
    obj_fun_norm = @(m_norm) compute_cost_lcurve(norm_to_true(m_norm), prior_params, xopt, yopt, zopt, insarx, insary, look, ...
        insar_lengths, nanstatbeginning, gps_obs, insar_obs, sigma, cinv_full, ...
        gps_weight, insar_weight, prior_weight, paramNames);
    
    % 4. Set up the [0, 1] bounds for the optimizer
    lb_norm = zeros(1, length(lb));
    ub_norm = ones(1, length(ub));
    
    % 5. Normalize the warm-start initial point
    x0_true = start_params;
    x0_norm = (x0_true - lb) ./ (ub - lb);
    
    % Safety catch: Ensure the starting point is strictly within [0, 1]
    x0_norm = max(0, min(1, x0_norm));
    
    % 6. Initialize Pattern Search
    options = optimoptions('patternsearch', ...
        'Display', 'off', ...           % Turned off so you don't get the confusing stop messages
        'UseParallel', true, ...
        'UseCompletePoll', true, ...
        'MeshTolerance', 1e-4, ...  
        'StepTolerance', 1e-4);
    
    % Run Pattern Search in the clean [0, 1] space
    [optParams_norm, ~] = patternsearch(obj_fun_norm, x0_norm, [], [], [], [], lb_norm, ub_norm, [], options);
    
    % Convert the optimal [0, 1] result back to physical parameters
    optParams_i = norm_to_true(optParams_norm);
    
    % 7. Re-evaluate at the optimal physical parameters to get isolated L2 metrics
    [~, gps_l2, insar_l2, prior_l2] = compute_cost_lcurve(optParams_i, prior_params, xopt, yopt, zopt, ...
        insarx, insary, look, insar_lengths, nanstatbeginning, gps_obs, insar_obs, ...
        sigma, cinv_full, gps_weight, insar_weight, prior_weight, paramNames);
end


function [cost, gps_l2, insar_l2, prior_l2] = compute_cost_lcurve(m, m_known, x, y, z, insarx, insary, look, insar_lengths, nanstatsbeg, ...
                                            gps_obs, insar_obs, sigma, cinv_full, gps_w, insar_w, prior_w, paramNames)
    % 1. Compute forward models exactly as the MCMC script does
    data_pred = create_MCMC_data(m, m_known, x, y, z, insarx, insary, look, insar_lengths, nanstatsbeg);
    
    gps_pred = data_pred(1:length(gps_obs));
    insar_pred = data_pred(length(gps_obs)+1:end);
    
    % 2. GPS L2 Metric (Sum of squared standardized residuals)
    gps_res = (gps_pred - gps_obs) ./ sigma;
    gps_l2 = sum(gps_res.^2);
    
    % 3. InSAR L2 Metric (Mahalanobis distance)
    insar_res = insar_pred - insar_obs;
    insar_l2 = insar_res' * cinv_full * insar_res;
    
    % 4. Prior L2 Calculation using your stored distributions
    persistent pDists
    if isempty(pDists)
        try
            temp = load('Data/paramDists.mat');
            pDists = temp.paramDists;
        catch
            pDists = struct();
        end
    end
    
    prior_l2 = 0;
    for idx = 1:length(paramNames)
        name = paramNames{idx};
        if isfield(pDists, name)
            dist = pDists.(name).dist;
            % Negative log-likelihood formulation
            p_val = pdf(dist, m(idx));
            if p_val > 0
                prior_l2 = prior_l2 - log(p_val);
            else
                prior_l2 = prior_l2 + 1e6; % large penalty to restrict search logic
            end
        end
    end
    L_scaling = (1/(length(gps_obs) + length(insar_obs)));
    % 5. Total composite cost
    cost = (gps_w * gps_l2 + insar_w * insar_l2) * L_scaling + prior_w * prior_l2;
end