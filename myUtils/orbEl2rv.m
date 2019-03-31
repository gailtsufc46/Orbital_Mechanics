function [r, v, h] = orbEl2rv(a, e, theta, OMEGA, omega, inc, mu)
% Written by Garrett Ailts
%
% Usage: [r, v, h] = orbEl2rv(a, e, theta, OMEGA, omega, inc)
%
% Description: Function takes a set of six keplerian orbital elements
% and computes the position and velocity vectors in the ECI frame. The
% specific angular momentum is also returned as an output.
%
% Inputs: a - semi major axis (km)
%         e - eccentricity
%     theta - true anomaly (rad)
%     OMEGA - right ascension of ascending node (rad)
%     omega - arguement of perigee (rad)
%       inc - inclination (rad)
%
% Outputs: r - position vector ECI frame (km)
%          v - velocity vector ECI frame (km/s)
%          h - specific angualr momentum (km^2/s)

%% Constants and Coordinates
pr = [cos(theta); sin(theta); 0];
pv = [-sin(theta); e+cos(theta); 0];
h = sqrt(mu*abs(a*(1-e^2)));

%% Compute r and v In the Perifocal Frame
rp = h^2/mu/(1+e*cos(theta))*pr;
vp = (mu/h)*pv;

%% Rotate r and v Into the ECI Frame
Q = angle2dcm(OMEGA,inc,omega,'ZXZ');
r = Q'*rp;
v = Q'*vp;

end
