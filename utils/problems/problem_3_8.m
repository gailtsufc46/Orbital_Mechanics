%% Problem 3.8
%
% A satellite is in earth orbit for which the perigee altitude is 200 km 
% and the apogee altitude is 600 km. Find the time interval during which 
% the satellite remains above an altitude of 400 km.
%
% Given:
%   Earth (mu=398600, Re=6378)
%   rP = Re+200, rA = Re+600
%   alt = 400 km
%
% Find:
%   The time interval where the satellite remains above 400 km alt.
%
%
% Method:
%   0. Compute a, e, T from rP, rA
%   1. Compute the true anomaly values for r = Re+400
%   2. Compute the corresponding eccentric anomaly values
%   3. Compute the corresponding mean anomaly values
%   4. Compute the elapsed time from the mean anomaly values

Re = 6378;
mu = 398600;

rP = Re+200;
rA = Re+600;

rC = Re+400;  % radial distance for altitude of 400 km

%%   0. Compute a, e from rP, rA
e = (rA-rP)/(rA+rP);
a = .5*(rA+rP);
T = sqrt(a^3/mu)*2*pi;

%%   1. Compute the true anomaly values for r = Re+400
cosTh = (a*(1-e^2)/rC - 1)/e;
th    = acos(cosTh);

%%   2. Compute the corresponding eccentric anomaly values
E = EccAnomFromTrueAnom( th, e );

%%   3. Compute the corresponding mean anomaly values
M = MeanAnomFromEccAnomE( E, e );

%%   4. Compute the elapsed time from the mean anomaly values
dT = T - 2*M/(2*pi)*T;

fprintf(1,'Elapsed time: %f min.\n',dT/60)

