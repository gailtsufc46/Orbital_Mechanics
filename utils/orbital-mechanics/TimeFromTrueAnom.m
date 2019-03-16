function t = TimeFromTrueAnom( e, tA, n )

% Compute the time elapsed since periapse crossing until the given true
% anomaly.
%
% USAGE:
%   t = TimeFromTrueAnom( e, tA, n );
%
% Inputs:
%   e     Eccentricity
%   tA    True anomaly (rad)
%   n     Mean orbit rarte (rad/s)
%
% Outputs:
%   t     Time elapsed since periapse crossing (sec)

E  = 2*atan( sqrt( (1-e)/(1+e) ) * tan(tA/2) );
Me = E - e*sin(E);
t  = Me/n;
