% Problem 3.10 from the textbook

%% GIVEN
%=====================================

% the orbit period in seconds
T = 14*3600; 

% perigee radius (km)
rp = 10000;

% times, measured from perigee crossing
t1 = 0.5*3600;  
t2 = 1.5*3600;

% Earth orbiting satellite
mu = 398600.44; % (km^3/s^2)

%=====================================

%% Compute constant orbit parameters  

% orbit rate
n = 2*pi/T;

% semi major axis
a = ( mu/n^2 )^(1/3);

% eccentricity
e = 1-rp/a;

%% Find true anomaly for Time 1 ...

% Mean anomaly for time 1
Me1 = n*t1;

% Eccentric anomaly for time 1
E1 = EccAnomFromMeanAnom( Me1, e, 1e-8 );

% True anomaly for time 1
th1 = 2*atan( sqrt((1+e)/(1-e)) * tan(E1/2) );

%% Find true anomaly for Time 2 ...

% Mean anomaly for time 2
Me2 = n*t2;

% Eccentric anomaly for time 2
E2 = EccAnomFromMeanAnom( Me2, e, 1e-8 );

% True anomaly for time 2
th2 = 2*atan( sqrt((1+e)/(1-e)) * tan(E2/2) );

%% Change in true anomaly
dTheta = th2 - th1;

%% Area swept during this time ...

% First compute angular momentum "h"
h = sqrt(a*mu*(1-e^2));

% rate of area sweep is constant, equal to half the angular momentum
dA = (t2-t1)*h/2;
