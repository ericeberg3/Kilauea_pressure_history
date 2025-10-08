% Define the objective function
objectiveFcn = @(x) (x.x1 - 3)^2 + (x.x2 - 2)^2 + (x.x3 - 5)^2;

% Define the variable ranges and initial guesses
variableRange = [
    optimizableVariable('x1', [-5, 10], 'Type', 'real', 'Transform', 'none');
    optimizableVariable('x2', [-5, 10], 'Type', 'real', 'Transform', 'none');
    optimizableVariable('x3', [-5, 10], 'Type', 'real', 'Transform', 'none')
];

% Perform Bayesian optimization
results = bayesopt(objectiveFcn, variableRange, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'MaxObjectiveEvaluations', 30, ...
    'IsObjectiveDeterministic', true);

% Display the best point found
bestPoint = results.XAtMinObjective;
disp('Best Point Found:');
disp(bestPoint);

% Plot the optimization results
figure;
plot(results);
title('Bayesian Optimization Results');
xlabel('Iteration');
ylabel('Objective Value');

% Extract the fitted Gaussian process model
gpModel = results.TrainedGP;

% Evaluate the posterior distribution at new points for each parameter
xNew1 = linspace(-5, 10, 100)';
xNew2 = linspace(-5, 10, 100)';
xNew3 = linspace(-5, 10, 100)';

[mu1, sigma1] = predict(gpModel, table(xNew1, ones(size(xNew1)), ones(size(xNew1)), 'VariableNames', {'x1', 'x2', 'x3'}));
[mu2, sigma2] = predict(gpModel, table(ones(size(xNew2)), xNew2, ones(size(xNew2)), 'VariableNames', {'x1', 'x2', 'x3'}));
[mu3, sigma3] = predict(gpModel, table(ones(size(xNew3)), ones(size(xNew3)), xNew3, 'VariableNames', {'x1', 'x2', 'x3'}));

% Plot the posterior mean and variance for each parameter
figure;
subplot(3,1,1);
plot(xNew1, mu1, 'b-', 'LineWidth', 2);
hold on;
fill([xNew1; flip(xNew1)], [mu1 + 2*sigma1; flip(mu1 - 2*sigma1)], 'b', ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none');
plot(xNew1, mu1 + 2*sigma1, 'r--', 'LineWidth', 1);
plot(xNew1, mu1 - 2*sigma1, 'r--', 'LineWidth', 1);
xlabel('x1');
ylabel('Objective Value');
title('Posterior Mean and Confidence Interval for x1');

subplot(3,1,2);
plot(xNew2, mu2, 'b-', 'LineWidth', 2);
hold on;
fill([xNew2; flip(xNew2)], [mu2 + 2*sigma2; flip(mu2 - 2*sigma2)], 'b', ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none');
plot(xNew2, mu2 + 2*sigma2, 'r--', 'LineWidth', 1);
plot(xNew2, mu2 - 2*sigma2, 'r--', 'LineWidth', 1);
xlabel('x2');
ylabel('Objective Value');
title('Posterior Mean and Confidence Interval for x2');

subplot(3,1,3);
plot(xNew3, mu3, 'b-', 'LineWidth', 2);
hold on;
fill([xNew3; flip(xNew3)], [mu3 + 2*sigma3; flip(mu3 - 2*sigma3)], 'b', ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none');
plot(xNew3, mu3 + 2*sigma3, 'r--', 'LineWidth', 1);
plot(xNew3, mu3 - 2*sigma3, 'r--', 'LineWidth', 1);
xlabel('x3');
ylabel('Objective Value');
title('Posterior Mean and Confidence Interval for x3');
