%% Problem 3.6
%
% If the eccentricity of the elliptical orbit is 0.3, calculate, 
% in terms of the period T, the time required to fly from P to B.
%
% Given:
%   e=0.3, theta=pi/2
%
% Find:
%   time to go from perigee to point B in terms of period T
%
%
% Method:
%   1. Compute the mean anomaly at point B, MB
%   2. The fraction of the orbit is MB/(2*pi) = t/T. Solve for t...

e = 0.3;
thB = pi/2;

%%   1. Compute the mean anomaly at point B, MB
EB  = EccAnomFromTrueAnom(thB,e);
MB  = MeanAnomFromEccAnomE(EB,e);

%%   2. The fraction of the orbit is MB/(2*pi) = t/T. Solve for t/T...

t_over_T = MB/2/pi


