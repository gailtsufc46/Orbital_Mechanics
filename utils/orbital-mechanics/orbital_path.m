function [r,x,y] = orbital_path( a, e, theta )
%
% Compute the orbital path in the perifocal frame given the true anomaly.
%
% Inputs: 
%   a         Semi major axis (km)
%   e         Eccentricity 
%   theta     True anomaly (rad)
%
% Outputs:
%   r         Radius (km)
%   x         x-coordinate in perifocal frame. +x axis towards periapsis
%   y         y-coordinate in perifocal frame. 
%



r = a*(1-e^2)./(1+e*cos(theta));
x = r.*cos(theta);
y = r.*sin(theta);
