function E = EccAnomFromTrueAnom( trueAnom, e )
%
% Compute the eccentric anomaly given the true anomaly
%
% Inputs: 
%   trueAnom  True anomaly (rad)
%   e         Eccentricity 
%
% Outputs:
%   E         Eccentric anomaly (rad)
%

E = 2*atan( sqrt((1-e)/(1+e)) * tan(trueAnom/2) );
