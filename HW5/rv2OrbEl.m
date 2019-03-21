function [a, e, theta, OMEGA, omega, inc, h] = rv2OrbEl(r,v,mu)
% Written by Garrett Ailts
% Usage: [a, e, theta, OMEGA, omega, inc, h] = rv2el(r,v)
% Description: Computes the set of Keplerian orbital elements for the orbit
% described by the a position and velocity vector in the ECI frame
% Inputs: r - position vector ECI frame (km)
%         v - velocity vector ECI frame (km/s)
%
% Outputs: a - semi major axis (km)
%          e - eccentricity
%      theta - true anomaly (rad)
%      OMEGA - right ascension of ascending node (rad)
%      omega - arguement of perigee (rad)
%        inc - inclination (rad)
%          h - specific angualr momentum (km^2/s)

%% Define Constants
rmag = norm(r);
vmag = norm(v);
vr = dot(v,r/rmag);
%% Define ECI Coordinate Vectors
I = [1 0 0]; K = [0 0 1]; % J = [0 1 0];
%% Calculate Ang. Momentum Vector, Eccentricity Vector, and Line of Nodes
hvec = cross(r,v);
N = cross(K,hvec);
evec = (1/mu)*((vmag^2-mu/rmag)*r-rmag*vr*v);
%% Calculate Ang. Momentum, Eccentricity, and Semi-Major Axis
h = norm(hvec);
e = norm(evec);
a = (h^2/mu)*(1/(1-e^2));
%% Calculate True Anomaly
theta = acos((1/e)*(h^2/rmag/mu-1));
%% Calculate RAAN, Arguement of Perigee, and inclination
OMEGA = dot(I,N)/norm(N);
omega = dot(N,evec)/norm(N)/e;
inc = dot(K,hvec)/h;
end
