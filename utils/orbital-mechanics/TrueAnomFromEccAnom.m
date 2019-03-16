function th = TrueAnomFromEccAnom( E, e )
%
% Compute the true anomaly from the eccentric anomaly, for eccentric orbits
%
% Inputs: 
%   E         Eccentric anomaly (rad)
%   e         Eccentricity (0 <= e < 1)
%
% Outputs:
%   th        True anomaly (rad)
%

th = 2*atan( sqrt((1+e)./(1-e)) .* tan(E/2) );

