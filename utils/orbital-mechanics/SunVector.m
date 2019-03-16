function u = SunVector( jD )

%   Compute the sun vector in the ECI frame. 
%
%   Form:
%   u = SunVector( jD )
%
%   References: The 1993 Astronomical Almanac

% Days from J2000.0
n = jD - 2451545.0;

% Mean anomaly
g = 357.528 + 0.9856003*n;

% Ecliptic longitude
lam = rem(280.460 + 0.9856474*n,360) + 1.915*sind(g) + 0.02*sind(2*g);

% Obliquity of ecliptic
obOfE = 23.439 - 4.00e-7*n;

% Equatorial rectangular coordinates of the Sun 
sLam = sind(lam);

u = [cosd(lam); cosd(obOfE).*sLam; sind(obOfE).*sLam];

