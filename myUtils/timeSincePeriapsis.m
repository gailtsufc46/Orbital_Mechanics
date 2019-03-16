function t = timeSincePeriapsis(n,e,theta)
% Written by Garrett Ailts
% Usage: t = timeSincePeriapsis(n,e,theta);
% Description: Computes the time elasped in seconds since periapsis for an 
% eccentric Keplerian orbit.
% Inputs: n - mean orbit rate (rad/s), e - eccentricity, theta - true 
% anomaly (rad)
% Outputs: t - time since periapsis (s)

%% Calculate Eccentric Anomaly
E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
%% Calculate Mean Anomaly
M = E-e*sin(E);
%% Calculate Time
t=M/n;
end