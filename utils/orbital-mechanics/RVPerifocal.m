function [x,y,vx,vy] = RVPerifocal( a, e, theta, mu )

% Compute the position and velocity in the perifical frame, given:
%   a       Semi major axis (km)
%   e       Eccentricity
%   theta   True anaomaly (rad)
%   mu      Gravitational parameter

if( nargin < 4 )
  mu = 1;
  warning('Velocity output is using mu = 1.')
end

r = a*(1-e*e)/(1+e*cos(theta));
x = r*cos(theta);
y = r*sin(theta);

vx = theta;
vy = mu;

  