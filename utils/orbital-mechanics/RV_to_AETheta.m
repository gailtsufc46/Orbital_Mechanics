function [a,e,theta,h] = RV_to_AETheta( r, v, mu )
%
% Given r and v vectors (position and velocity), compute a, e, theta.
%
% Inputs: 
%   r     (3,1)   Position vector, km
%   v     (3,1)   Velocity vector, km
%   mu    (1,1)   Gravitational parameter, km^3/s^2
%
% Outputs:
%   a     (1,1)   Semi major axis, km
%   e     (1,1)   Eccentricity
%   theta (1,1)   True anomaly, rad
%

% use Earth's gravitational constant if none is provided
if( nargin < 3 )
  mu = 398600.44;
end
  
rvec = r;
vvec = v;

% vector magnitudes
r = sqrt(r'*r);
v = sqrt(v'*v);

% radial velocity
vr = vvec'*rvec/r;

% angular momentum
hvec = cross(rvec,vvec);
h = sqrt(hvec'*hvec);

% true anomaly
theta = atan2( h*vr*r, h^2-mu*r );

% eccentricity
e = h*vr/(mu*sin(theta));

% semi major axis
a = h^2/mu / (1-e^2);


