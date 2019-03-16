function theta = TrueAnomFromHypEccAnom( F, e )
%
% Compute the true anomaly given the hyperbolic eccentric anomaly 
%
% Inputs: 
%   F         Hyperbolic eccentric anomaly (rad)
%   e         Eccentricity (e>1)
%
% Outputs:
%   theta     True anomaly (rad)
%

theta = 2*atan( sqrt( (e+1)/(e-1) )* tanh(F/2) );
