function  [minit] = initialize_walkers(W, M, x0, xbnds, data_vec, data_sigma)
% This function initializes walkers by first randomly create 5* W walkers
% centered around x0. Then we pick out W walkers with the best liklihoods
% and use them as the initial ensemble of workers. 


% First we need to initialize more than enough walkers to pick out W
% walkers' initial poisitions. Here to be safe we initialize 5 * W walkers

minit_prop = zeros(M, 5*W); % preallocate memory for minit_prop

for w = 1:5*W
    for m = 1:M
        distance =  (-1+2*rand(1))*((xbnds(m,2)-xbnds(m,1))/10); % We want the initial walkers to be a tight ball around a good fitting set of parameters
        minit_prop(m,w) = x0(m) + distance; % a Mx5*W matrix of initial values 
    end
end

% Then we do 1 iteration with 4*W walkers. This step is not done in parallel
itr = 3* 5* W;
logPfuns = {@(m)logprior(m, xbnds) @(m)loglike(m, data_vec, data_sigma)};
[models,logP] = gwmcmc(minit_prop,logPfuns, itr, 'StepSize', 2.5, 'ThinChain', 1, 'Parallel', true, 'FileName', 'initial_models_logP', 'OutputData', 1, 'NCores', 4); 
logP = logP(1,:,:) + logP(2,:,:); % we want probability = log(likihood) + log(prior)
[maxk_like, loc_maxk] = maxk(logP(:,:,3), W); % loc_maxk marks the indices of walkers at itr = 3 with the best likihood
I = find(maxk_like > -1e5); % make sure that all the chosen likihoods are reasonable
loc_maxk = loc_maxk(I);
disp(['Expected ' num2str(W), ' walkers. Initiated with ', num2str(length(loc_maxk)), ' walkers'])

% Then we identify the models associated with walkers with high liklihoods
minit = models(:, loc_maxk, 3); % a MxW matrix of initial values 

save('initial_walkers.mat', 'minit')

delete(gcp('nocreate')) % terminate the existing parallel session
end