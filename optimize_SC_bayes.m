function optParams = optimize_SC_bayes(mSCguess, m_known, lb, ub, xopt, xtilt, yopt, ytilt, zopt, u1d, ...
    invStdPWRL, tiltstd, tiltreduced, nanstatbeginning)

priormeans = (ub + lb)/2;
priorvariances = (ub - lb)/3;

f = @(m)green_residuals(m, m_known, [xopt, xtilt], [yopt, ytilt], [zopt, 0], u1d, ...
    [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3),1/((tiltstd)), 1/((tiltstd))], tiltreduced(1:2), nanstatbeginning, priormeans, priorvariances);


params = [
    optimizableVariable('dpHMM', [lb(1), ub(1)], 'Type', 'real', 'Transform', 'none');
    optimizableVariable('vert_sd', [lb(2), ub(2)], 'Type', 'real', 'Transform', 'none');
    optimizableVariable('horiz_sd', [lb(3), ub(3)], 'Type', 'real', 'Transform', 'none');
    optimizableVariable('dip', [lb(4), ub(4)], 'Type', 'real', 'Transform', 'none');
    optimizableVariable('dpSC', [lb(5), ub(5)], 'Type', 'real', 'Transform', 'none');
];

% Perform the optimization
results = bayesopt(@(m) f([m.dpHMM, m.vert_sd, m.horiz_sd, m.dip, m.dpSC]), ...
                   params, ...
                   'MaxObjectiveEvaluations', 200, 'UseParallel',false, 'GPActiveSetSize', 400, 'IsObjectiveDeterministic', false);


optParams = results.XAtMinObjective;
disp('Optimal Parameters:');
disp(optParams);