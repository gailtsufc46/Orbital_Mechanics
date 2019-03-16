%% Problem 6.29 from the textbook
% Apse line rotation
%
% Note: first part similar to problem 3.18

%% constants
mu = 398600.44; % Earth gravitational parameter (km^3/s^2)
Re = 6378.14;   % Earth equatorial radius (km)

%% given
rA = 12756;           % initial radius
vA = 6.5992;          % initial velocity (total)
gammaA = 20*pi/180;   % initial flight path angle
thB = 150*pi/180;     % true anomaly where impulsive maneuver occurs
dVT = 0.7582;         % Transverse velocity change at B
dVR = 0;              % Radial velocity change at B

%% Initial computations...

% transverse velocity and radial velocity
vtrans = vA*cos(gammaA);
vrad = vA*sin(gammaA);

% angular momentum
h = rA*vtrans;

% orbit equation: r = (h^2/mu) / (1 + e*cos(th))
% rearrange to solve for e*cos(th)
ecosth = h^2 / mu /rA - 1;

% radial velocity equation: vrad = (mu/h)*e*sin(th)
% rearrange to solve for e*sin(th)
esinth = vrad*h/mu;

%% solve for th0, e, E, M, t

% tan(th) = sin(th) / cos(th) = [e*sin(th)] / [e*cos(th)]
% solve for true anomaly 
% let's try it THE CORRECT way now...
thA = atan2(esinth,ecosth);

% solve for eccentricity
e = esinth/sin(thA);

% semi major axis 
a = (h^2/mu)/(1-e^2);

%% Time of flight from A to B
tA = TimeSincePeriapsis(thA,a,e,mu);
tB = TimeSincePeriapsis(thB,a,e,mu);

tOF = (tB-tA)/3600

%% Solve for new true anomaly at B after apse line rotation maneuver
rB      = h^2/mu/(1+e*cos(thB));
vRadB   = mu/h*e*sin(thB);
vTransB = mu/h*(1+e*cos(thB));
tanThB2 = vTransB^2/(mu/rB) * ...
  (vTransB+dVT)*(vRadB+dVR) / ...
  ( (vTransB+dVT )^2*e*cos(thB) + (2*vTransB+dVT)*dVT );

% two possible solutions:

if( tanThB2<0 )
  thB2_neg = atan(tanThB2);
else
  thB2_neg = atan(tanThB2)-pi;
end
thB2_pos = thB2_neg+pi;

% determine which true anomaly is correct based on sign of radial velocity
if( vRadB > 0 )
  thB2 = thB2_pos;
else
  thB2 = thB2_neg;
end

%% Apse Line Rotation
eta = (thB2-thB)*180/pi
if( eta<0 )
  disp('Counter-Clockwise')
else
  disp('Clockwise')
end



