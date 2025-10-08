function t=  travel_time(hypo, sta, vel);
% t=  travel_time(hypo, sta, vel);
%
% compute travel time given hypocentral parameters and station locations,
% assuming straight ray paths, and uniform half-space with constant
% velocity, vel
%
% Imput:
%       hypo = hypocentral coords(4,N), [East, North, Height, origin_time]
%       sta  = station coordinates (3,M) [East, North, Height]
%       vel  = wave speed
%
%   Note source depth is a negative height.

Nsta = size(sta,2);
vecs = repmat(hypo(1:3,:), 1, Nsta) - sta;
dist = sqrt(vecs(1,:).^2 + vecs(2,:).^2 + vecs(3,:).^2);
t    = dist/vel + hypo(4,:);

