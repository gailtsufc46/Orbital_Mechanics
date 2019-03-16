%% Problem 3.18 from the textbook

%% Problem statement
% An incoming object is sighted at a radius of 100,000 km with a speed
% of 6 km/s and a flight path angle of 80°. 
% (a) Will it impact the earth or fly by? 
% (b) What is the time either to impact or closest approach?

%% constants
mu = 398600.44; % Earth gravitational parameter (km^3/s^2)
Re = 6378.14;   % Earth equatorial radius (km)

%% given
r0 = 100000;          % initial radius
v0 = 6;               % initial velocity (total)
gamma0 = -80*pi/180;  % initial flight path angle

%% Initial computations...

% transverse velocity and radial velocity
vtrans = v0*cos(gamma0);
vrad = v0*sin(gamma0);

% angular momentum
h = r0*vtrans;

% orbit equation: r = (h^2/mu) / (1 + e*cos(th))
% rearrange to solve for e*cos(th)
ecosth = h^2 / mu /r0 - 1;

% radial velocity equation: vrad = (mu/h)*e*sin(th)
% rearrange to solve for e*sin(th)
esinth = vrad*h/mu;

%% method 1 ... This method has a mistake in it
% solve for th0, e, F, Mh, t

% tan(th) = sin(th) / cos(th) = [e*sin(th)] / [e*cos(th)]
% solve for true anomaly 
% let's try it this way first...
th0 = atan(esinth/ecosth);
% the spacecraft is approaching Earth, so th0 should be <0
if( th0>0 )
  th0 = -th0; % This may seem okay, but its actually wrong!
end

% solve for eccentricity
e = esinth/sin(th0);

% semi major axis 
a = (h^2/mu)/(e^2-1);

% what is the asymptote angle?
thinf = acos(-1/e);

% solve for hyperbolic eccentric anomaly, F
F = HypEccAnomFromTrueAnom(th0,e);

% solve for hyperbolic mean anomaly, Mh
Mh = e*sinh(F)-F;

% solve for "time since perigee crossing"
t = Mh / ( h^3 / mu^2 / (e^2-1)^1.5 );

% the time in hours
tHr = t/3600;

% Check this result... 
% Use common sense. Is it feasible? Why or why not?
% Use the orbit equation with th0 to check the initial radius
r_check1 = h^2/mu / (1+e*cos(th0)); % should be same as r0


%% method 2 ... This method is correct
% solve for th0, e, F, Mh, t

% tan(th) = sin(th) / cos(th) = [e*sin(th)] / [e*cos(th)]
% solve for true anomaly 
% let's try it THE CORRECT way now...
th02 = atan2(esinth,ecosth);

% solve for eccentricity
e2 = esinth/sin(th02);

% what is the asymptote angle?
thinf2 = acos(-1/e2);

% solve for hyperbolic eccentric anomaly, F
F2 = HypEccAnomFromTrueAnom(th02,e2);

% solve for hyperbolic mean anomaly, Mh
Mh2 = e*sinh(F2)-F2;

% solve for "time since perigee crossing"
t2 = Mh2 * h^3 / mu^2 / (e2^2-1)^1.5;

% the time in hours
tHr2 = t2/3600;

% Check this result... 
% Use the orbit equation with th0 to check the initial radius
r_check2 = h^2/mu / (1+e*cos(th02)); % should be same as r0

