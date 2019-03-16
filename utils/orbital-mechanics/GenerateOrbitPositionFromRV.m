function [xI,yI,zI] = GenerateOrbitPositionFromRV( r0, v0, mu )

% Given the position and velocity vectors at some point in time (r0,v0),
% and the gravitational parameter of the central body (mu), compute the
% x,y,z position data for the Keplerian orbit. 
%
% If the orbit is elliptic or circular (e<1), we compute position data over
% the full range of true anomaly angles, from 0 to 2*PI.
%
% If the orbit is parabolic or hyperbolic (e>=1), then we compute position
% data over a feasible range of true anomaly angles, from:
%   -thetaLim+TOL to +thetaLim-TOL
%
% where "thetaLim" is the true anomly limit of acos(-1/e), and TOL>0 is an
% arbitrary small tolerance in radians.

[a,i,W,w,e,th0] = OrbitalElementsFromRV(r0,v0,mu);

if( e<1 )
  th = linspace(0,2*pi,300);
else
  thLim = acos(-1/e);
  th = linspace(-thLim+.1, thLim-.1, 300 );
end

% first compute position in perifocal coordinates
r = a*(1-e^2)./(1+e*cos(th));
x = r.*cos(th);
y = r.*sin(th);

if( isnan(W) )
  W = 0;
end

if( isnan(w) )
  w = 0;
end

% next rotate to ECI coordinates...

% this rotates from ECI to perifocal
Q = [cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1] * ...
  [1 0 0;0 cos(i) sin(i); 0 -sin(i) cos(i)] * ...
  [cos(W) sin(W) 0; -sin(W) cos(W) 0; 0 0 1];

% the transpose rotates from perifocal to ECI
M = Q';              

rECI = M*[x;y;zeros(size(x))];
xI = rECI(1,:);
yI = rECI(2,:);
zI = rECI(3,:);
