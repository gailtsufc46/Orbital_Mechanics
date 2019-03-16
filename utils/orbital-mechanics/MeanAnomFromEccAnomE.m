function Me = MeanAnomFromEccAnomE( E, e )
%
% Compute mean anomaly for a given eccentric anomaly in an eccentric orbit
%
% Inputs: 
%   E         Eccentric anomaly
%   e         Eccentricity (0 <= e < 1)
%
% Outputs:
%   Me        Mean anomaly
%

Me = E - e*sin(E);
