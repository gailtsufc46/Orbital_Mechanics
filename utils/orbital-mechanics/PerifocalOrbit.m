function [r,x,y] = PerifocalOrbit( a, e, theta )

%PerifocalOrbit  Orbit trajeectory in perifocal frame.
%
% Compute the orbital path in the perifocal frame given the true anomaly.
%
% [r,x,y] = PerifocalOrbit( a, e, theta )
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
%
%   Copyright 2014 Joseph Mueller


r = a*(1-e^2)./(1+e*cos(theta));
x = r.*cos(theta);
y = r.*sin(theta);
