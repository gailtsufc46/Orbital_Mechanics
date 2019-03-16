function t = TimeSincePeriapsis( trueAnom, a, e, mu )
%
% Compute the time elapsed from periapsis crossing to given true anomaly
%
% Inputs: 
%   trueAnom  True anomaly (rad)
%   a         Semi major axis 
%   e         Eccentricity 
%   mu        Gravitational parameter (km^3/s^2) 
%
% Outputs:
%   t         Time measured since periapsis crossing (seconds) 
%

E = EccAnomFromTrueAnom( trueAnom, e );
Me = MeanAnomFromEccAnomE( E, e );
n = OrbRate( a, mu );
t = Me/n;
