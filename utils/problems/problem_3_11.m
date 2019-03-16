%% Problem 3.11

% A satellite in earth orbit has perigee and apogee radii of:
%   rp = 7500 km 
%   ra = 16,000km
% 
% Find its true anomaly 40 minutes after passing true anomaly of 80°.

% Method:
% 0. Compute a, e, n for the orbit
% 1. Compute E at theta = 80
% 2. Compute M at theta = 80
% 3. Compute delta M for 40 minutes
% 4. Compute M at new time (40 minutes after theta = 80)
% 5. Compute E at new time (40 minutes after theta = 80)
% 6. Compute theta at new time (40 minutes after theta = 80)

%% Given
rp = 7500;
ra = 16000;
mu = 398600.;

%% 0. Compute a, e, n for the orbit
e = (ra-rp)/(rp+ra);
a = .5*(rp+ra);
n = sqrt(mu/a^3);


%% 1. Compute E at theta = 80
th80  = 80*pi/180;
E80   = EccAnomFromTrueAnom(th80,e)

%% 2. Compute M at theta = 80
M80   = MeanAnomFromEccAnomE(E80,e)

%% 3. Compute delta M for 40 minutes
dT    = 40*60;  % elapsed time in seconds
dM    = dT*n   % change in mean anomaly over this time

%% 4. Compute M at new time (40 minutes after theta = 80)
M2    = M80 + dM

%% 5. Compute E at new time (40 minutes after theta = 80)
E2    = EccAnomFromMeanAnom(M2,e,1e-8)

%% 6. Compute theta at new time (40 minutes after theta = 80)
theta2 = TrueAnomFromEccAnom(E2,e)

fprintf(1,'True anomaly 40 minutes later is: %f\n',theta2*180/pi);
