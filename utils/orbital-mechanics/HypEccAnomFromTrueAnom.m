function F = HypEccAnomFromTrueAnom( theta, e )
%
% Compute the hyperbolic eccentric anomaly given the true anomaly
%
% Inputs: 
%   theta     True anomaly (rad)
%   e         Eccentricity 
%
% Outputs:
%   F         Hyperbolic eccentric anomaly (rad)
%

F = 2*atanh( sqrt( (e-1)/(e+1) )* tan(theta/2) );
