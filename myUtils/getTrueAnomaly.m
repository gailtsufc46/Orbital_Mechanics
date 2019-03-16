function nu = getTrueAnomaly(n,e,t)
% Written by Garrett Ailts
% Usage: nu = getTrueAnomaly(n,e,t); 
% Description: Calculate the true anomaly of an eccentric Keplerian orbit
% Inputs: n - mean orbit rate (rad/s), e - eccentricity, t - time since
% periapsis (s)
% Ouputs: nu - true anomaly in radians


%% Calculate Mean Anomaly
M = n*t;
%% Calculate Eccentric Anomaly
E = mean2Ecc(e,M);
%% Calculate True Anomaly
nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end
