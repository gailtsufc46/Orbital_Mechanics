function T = OrbPeriod( a, mu )
%
% Compute the orbit period given the semi major axis and grav. param. mu
%
% Inputs: 
%   a         Semi major axis 
%   mu        Gravitational parameter (km^3/s^2) 
%
% Outputs:
%   T         Orbital period (sec)
%

T = 2*pi*sqrt(a^3 / mu );
