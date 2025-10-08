%% This function estimates RW noise using a quiet period
function rw_stddev = GetRandomWalk(sdh_clean)
    rw_data = sdh_clean.d(17818:52599, 1);
    
    % Sample rate of tiltmeter is 1/min, so sqrt(delta t) = sqrt(1 min)
    % Thus, the random walk std is simply the sqrt of the differenced variance
    rw_stddev = std(diff(rw_data)) * sqrt(60); % must mult. by sqrt(60) to convert from sqrt(min) -> sqrt(hr)
    disp("RW Standard Deviation = " + rw_stddev + " From " + datestr(sdh_clean.t(5.24e4)) + " To " + datestr(sdh_clean.t(5.38e4)))
    
    %% Plotting frequency spectrum
    
    fs = 1/((sdh_clean.t(2) - sdh_clean.t(1)) * 86400); % sampling freq (Hz)
    
    [pxx,f] = pwelch(rw_data, [], [], [], fs); 
    
    % We want to fit between fLow = 1e-5 Hz and fHigh = 1e-2 Hz
    idxFit = f >= 1e-5 & f <= 1e-2;
    freqToFit = f(idxFit);
    psdToFit  = pxx(idxFit);
    
    % Define the model as a function handle
    modelFun = @(A, f) A ./ (f.^2);
    
    % Initial guess for A
    A0 = 1;
    
    % Use lsqcurvefit
    options = optimset('Display','off');  % optional
    A_fit = lsqcurvefit(modelFun, A0, freqToFit, psdToFit, [], [], options);
    
    X = log(freqToFit);
    Y = log(psdToFit);
    
    p = polyfit(X, Y, 1);   % Linear fit
    slope = p(1);
    intercept = p(2);
    
    % slope should be close to -2 if it's truly 1/f^2
    % intercept = ln(A) => A = exp(intercept)
    A_fit_log = exp(intercept);
    
    X = log(freqToFit);
    Y = log(psdToFit);
    
    p = polyfit(X, Y, 1);   % 1 => linear fit
    slope = p(1);
    intercept = p(2);
    
    % slope should be close to -2 if it's truly 1/f^2
    % intercept = ln(A) => A = exp(intercept)
    A_fit_log = exp(intercept);
    
    % Plotting code:
    % figure(11);
    % loglog(freqToFit, psdToFit, 'bo', 'MarkerSize', 5);
    % hold on;
    % 
    % % If you used lsqcurvefit (linear domain fitting):
    % psdFitted = A_fit ./ (freqToFit.^2);
    % loglog(freqToFit, psdFitted, 'r-', 'LineWidth', 2);
    % 
    % % (Or if you used the polyfit approach)
    % % psdFitted = A_fit_log ./ (freqToFit.^2); % slope ~ -2
    % % loglog(freqToFit, psdFitted, 'r-', 'LineWidth', 2);
    % 
    % grid on;
    % xlabel('Frequency (Hz)');
    % ylabel('PSD');
    % legend('Data', '1/f^2 Fit');
    % title('1/f^2 Fit to PSD');
end