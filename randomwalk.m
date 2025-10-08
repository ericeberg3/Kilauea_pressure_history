downsamplerate = 720;
load(['data/u_mm_' int2str(downsamplerate) '.mat']);
load('data/gps_time_series.mat');
load('data/GPS_seism_locations.mat', 'GPSNameList', 'GPS_llh');
load('data/hawaii_line_new.mat', 'coast_new', 'new_pit')
load('data/sdh_clean.mat')

t = decyear(sdh_clean.t);
year = floor(t);
partialYear = mod(t,1);
date0 = datenum(num2str(year),'yyyy');
date1 = datenum(num2str(year+1),'yyyy');
daysInYear = date1 - date0;
t = date0 + partialYear .* daysInYear;
t = datetime(t, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy');
% collapset = datetime(collapset, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy');
sigma = var(diff(sdh_clean.d(1e4:6e4, 1)));
rwsigma = sigma * sqrt(find(t == datetime('2018-08-10')) - find(t == datetime('2018-06-09')));
% disp(rwsigma)

% S = sparse(tril(ones(5e4-1, 5e4-1)));
% C = cov(sdh_clean.d(1e4:6e4, 1));
lambda = 1;
% sigmasq = S \ S' \ C;
% disp(sigmasq)

figure(1);
% Use feb 08 to Mar 15 for random walk fitting
clf;
hold on;
plot(t(1e4:6e4), sdh_clean.d(1e4:6e4, 1));


%% Generate a random walk using cumsum
% % Example data
distances = ones(1, 100);  % distances between benchmarks (in meters)
sigma = 1e-2;  % measurement precision (standard deviation)
height_differences = normrnd(0,sigma, 1, 5e4);  % observed height differences (in meters)

% Number of benchmarks
n = length(distances) + 1;

% Initialize the covariance matrix for section height differences
C_delta_h = diag(sigma^2 * distances);

% Construct the lower triangular matrix S
S = tril(ones(n-1, n-1));

% Calculate the covariance matrix for benchmark elevations
C_H = S * C_delta_h * S';

% Calculate the benchmark elevations by summing the section height differences
H = cumsum([0, height_differences]) - 376.904;

plot(t(1e4:6e4), H);
