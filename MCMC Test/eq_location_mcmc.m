% choose some stations
clearvars
rng('default');  % for testing


Nsta = 8;
sta  = zeros(3,Nsta);
sta(1:2,:) = 2*(rand(2,Nsta)-0.5);  % locations from -1 to 1

% try making deep stations
%sta(3,:) = -rand(1,Nsta);

% make true hypocenter
hyp_true = [0; 0; -1; 0];
xbnds  = [-1 1;-1 1;-3 0 ;-1 1];

%synthetic data
vel = 1;
t =  travel_time(hyp_true, sta, vel);
t_true = t;

% add noise
sigma = 0.01*mean(t);  % for now known error
t = t + sigma*randn(size(t));


%% to run mcmc
x0 = [0.3; 0.3; -0.8; 0.1];
xstep = 0.05;   % step size  step size in t
Niter=10000;



[x_keep, L_keep, count] = mcmc('travel_time',t,x0,xstep,xbnds,sigma,Niter,sta,vel);

%% plot some of the output

figure
subplot(311)
%plot(log(L_keep)); box on
plot(L_keep); box on
%ylim([-50 50])
ylabel('Log(P)')
subplot(312)
plot(x_keep(1,:));
grid on; box on
ylim([-1 1]), ylabel('X-coord')
subplot(313)
plot(x_keep(3,:)); 
grid on; box on
xlabel('Iterations'); ylabel('Z-coord'); %ylim([-1 1])
accratio = ['Acceptance rate = ' num2str(count/Niter)];
subplot(311)
title(accratio)


%% plot results
burn = 1000;  % need to find by eye

%MAP estimate
[~,I] = max(L_keep);
x_keep(:,I);
map_hyp = x_keep(:,I);
tpred =  travel_time(x_keep(:,I), sta, vel);
r = t - tpred;

figure
%plot3(x_keep(1,:), x_keep(2,:), x_keep(3,:),'k.');
hold on

plot3(x_keep(1,burn:end), x_keep(2,burn:end), x_keep(3,burn:end), '.b');
plot3(x_keep(1,1:burn), x_keep(2,1:burn), x_keep(3,1:burn), '.k');
plot3(hyp_true(1), hyp_true(2), hyp_true(3), 'r+', 'LineWidth',5);
plot3(map_hyp(1), map_hyp(2), map_hyp(3), 'r+', 'LineWidth',5);

grid on; box on;
xlim([-1 1]); ylim([-1 1]); zlim([-2 0])

% plot stations
plot3(sta(1,:), sta(2,:), sta(3,:), 'r^')

%% plot historgrams
figure
subplot(221); hist(x_keep(1,burn:end)), xlabel('x-coord')
subplot(222); hist(x_keep(3,burn:end)), xlabel('z-coord')
subplot(223); hist(x_keep(4,burn:end)), xlabel('origin t')
subplot(224); plot(x_keep(3,burn:end), x_keep(4,burn:end), '.')
xlabel('z-coord'); ylabel('origin t')

%MAP estimate
[~,I] = max(L_keep);
x_keep(:,I)
tpred =  travel_time(x_keep(:,I), sta, vel);
r = t - tpred
return

title('Event location');
xlabel('Coordinate x (km)'); ylabel('Coordinate y (km)');
legend('MCMC', 'No Burn In', 'True solution');

