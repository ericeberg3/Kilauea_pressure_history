function  [x_keep, L_keep, count, gps_l2, insar_l2, prior_l2] = mcmc(func,data,x0,xstep,xbnds,sigma,cinv_full,Niter, ...
    gps_weighting,prior_weight,priorNames,burn,varargin)
%
% [x_keep, L_keep, count] = mcmc(func,data,x0,xstep,sigma,Niter,varargin)
%
% subroutine for MCMC sampling using Metropolis-Hasting w/ normal
% distribution.
%
% Inputs:
%       func = function name eg:   mcmc('travel_time',....)
%               parameter required for function in varargin
%       data = vector of observations
%       x0   = initial estimate of parameter vector
%       xstep = step size in all parameter directions
%       xbnds = bounds (
%       sigma = sigma of normal distribution 
%       Niter = number of iterations
%
% Outputs:
%       x_keep = array of samples
%       L_keep = likelihood of samples
%       count  = number of accepted. Acceptance ratio is count/Niter
%
%  P Segall: 2012

fun = fcnchk(func);
load Data/paramDists.mat paramDists;
% fun = str2func(func);

%number of elements in x
Nparams = length(x0);
% check dimensions of bounds
if( size(xbnds,1) ~= Nparams |  size(xbnds,2) ~= 2)
    disp('Dimension of xbnds is not valid')
    return
end

x_keep=zeros(Nparams,Niter); 
L_keep=zeros(1,Niter); 
N_gps = length(varargin{3})*3;
N_insar = length(varargin{6});

x = x0;
% prior_weight = 0*5e2;
dprop = fun(x, varargin{:});
L_scaling = (1/length(data)) * 0.5;

L_gps = -0.5 * sum(((data(1:N_gps) - dprop(1:N_gps))./sigma(1:N_gps)).^2);
L_insar = -0.5 * (data(N_gps+1:end) - dprop(N_gps+1:end))' * cinv_full * (data(N_gps+1:end) - dprop(N_gps+1:end));
prior_prob = 0;
for j = 1:numel(priorNames)
    name = priorNames{j};
    prior_prob = prior_prob + log(pdf(paramDists.(name).dist, x(j)));
end
L = L_scaling * (L_gps*gps_weighting + L_insar + prior_weight*prior_prob);
barLength = 25;
prevStr = '';

count=0;
xstep_int = xstep;
xstep = max(abs(xbnds(:,2)'), abs(xbnds(:,1)')) .* xstep_int;

p          = numel(x0);
Sraw       = 1e-4 * eye(p);     % running (un‑scaled) covariance
log_s      = log(0.1);          % log‑scale parameter  (starts small)
targetAcc  = 0.234;             % optimal random‑walk rate
Y          = 0.05;              % adaptation gain (decays)
eps_reg    = 1e-10;             % jitter

mean_prop  = x0;
accCount   = 0;

% Set step for negative values to be lower bounds:
% xstep(1) = xbnds(1,1) * xstep_int(1); % dpHMM
% xstep(4) = xbnds(4,1) * xstep_int(4); % dpSC

for k=1:Niter
    % generate proposal
    % xprop = x + xstep.*(rand(Nparams,1)'-0.5);

    % scale   = exp(log_s);
    % SigProp = Sraw * scale^2 + eps_reg*eye(p);
    xprop = x + xstep.*(rand(Nparams,1)'-0.5);

    % HMM_volume = xprop(2);
    % opt_vert_sd = (3/(4*pi) * HMM_volume * (aspect_ratio_HMM^2))^(1/3);
    % surf_depth = xprop(5) - opt_vert_sd;
    % check bounds
    if all(xprop > xbnds(:,1)') && all(xprop < xbnds(:,2)')
        % && xprop(2) < xprop(1) && xprop(end) < xprop(end-1) % GPS pressure < insar pressure
        
        % prediction
        dprop = real(fun(xprop, varargin{:}));
        
        % likelihood
        % Lprop = prod(normpdf( data-dprop, 0, sigma));
        %Lprop = exp(-0.5*(norm(data-dprop))^2/sigma^2);
        % log likelihood
        
        % Maybe make this L1 norm instead of L2
        % Look at Kyle's distributions and see how smooth they are
        % Try MC hammer algorithm to parallelize 
        % Should be running 1M simulations
        prior_prob = 0;
        for j = 1:numel(priorNames)
            name = priorNames{j};
            prior_prob = prior_prob + log(pdf(paramDists.(name).dist, xprop(j)));
        end


        L_gps = -0.5 * sum(((data(1:N_gps) - dprop(1:N_gps))./sigma(1:N_gps)).^2);
        L_insar = -0.5 * (data(N_gps+1:end) - dprop(N_gps+1:end))' * cinv_full * (data(N_gps+1:end) - dprop(N_gps+1:end));

        Lprop = L_scaling * (L_gps*gps_weighting + L_insar + prior_weight*prior_prob);

        u=rand(1,1);

        %if (L==0 || u <= Lprop/L)
        if (log(u) <= Lprop - L)
             count=count+1;
             x = xprop;
             L = Lprop;
        end 
        accCount = accCount + 1;
    end

    x_keep(:,k) = x;
    L_keep(k) = L;

    % if k > 500          % first 500 steps = burn‑in, no adaptation
    %     alpha = 1/(k-500);  % Robbins‑Monro weight
    %     % 3a update mean & raw covariance
    %     oldMean = mean_prop;
    %     mean_prop = mean_prop + alpha*(x - mean_prop);
    %     Sraw      = (1-alpha)*Sraw + alpha*((x - oldMean)'*(x - oldMean));
    % 
    %     % 3b update global scale toward target acceptance
    %     log_s = log_s + Y*alpha*( (accCount/k) - targetAcc );
    % end

    if(mod(k, 0.2e3) == 0)
        % --- Update Progress Bar ---
        percent = k / Niter;
        acc_ratio = count/k;
        numEquals = round(percent * barLength);
        bar = [repmat('=', 1, numEquals) repmat(' ', 1, barLength - numEquals)];
        progressStr = sprintf('[%s] %3.1f%%', bar, acc_ratio*100);
        
        % Erase the previous progress bar
        if ~isempty(prevStr)
            fprintf(repmat('\b', 1, length(prevStr)));
        end
        
        % Print the new progress bar and store the string for the next iteration
        fprintf('%s', progressStr);
        prevStr = progressStr;
        
        drawnow;  % force the update
    end
end
fprintf('\n');
disp("Acceptance ratio: " + count/Niter)

% p        = 0.01;                              % 5 % neighbourhood
% thresh   = prctile(real(L_keep(burn:end)), 100*(1-p));
% idxTop   = L_keep(burn:end) >= thresh;
% x_keep_burnt = x_keep(:, burn:end);
% optParams = mean(real(x_keep_burnt(:, idxTop)), 2);

[~, ind] = max(L_keep(burn:end));
x_keep_burnt = x_keep(:, burn:end);
optParams = x_keep_burnt(:, ind);

dprop = fun(optParams', varargin{:});

gps_l2 = sum(((data(1:N_gps) - dprop(1:N_gps))./sigma(1:N_gps)).^2);
insar_l2 = (data(N_gps+1:end) - dprop(N_gps+1:end))' * cinv_full * (data(N_gps+1:end) - dprop(N_gps+1:end));
prior_l2 = 0;
for j = 1:numel(priorNames)
    name = priorNames{j};
    prior_l2 = prior_l2 + log(pdf(paramDists.(name).dist, optParams(j)));
end

end