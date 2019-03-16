%% Problem 3.14

% Calculate the time required for a spacecraft launched into a parabolic 
% trajectory at a perigee altitude of 200 km to leave the earth?s sphere of 
% influence (see Table A.2).
% 
% Sphere of influence for Earth is rSOI = 925000 km.

% Given:
%   Earth (mu=398600, Re=6378)
%   Parabola (e=1)
%   Perigee radial distance: rp = Re+200
%   Sphere of influence radius: rSOI = 925000
%
% Method:
% 0. Compute velocity at perigee knowing orbit energy is zero
% 1. Compute angular momentum h at perigee with rp, vp
% 2. Compute true anomaly at rSOI using orbit equation
% 3. Solve for parabolic mean anomaly using the cubic expression for it
% 4. Compute the elpased time from the parabolic mean anomaly

%% Given
mu = 398600.;
Re = 6378.;
e = 1.;
rp = Re+200;
rSOI = 925000;

%% 0. Compute velocity at perigee knowing orbit energy is zero
vp = sqrt( 2*mu/rp );

%% 1. Compute angular momentum h at perigee with r, v
h = rp*vp;

%% 2. Compute true anomaly at rSOI using orbit equation
cosTh = (h^2/mu)/rSOI - 1;
theta = acos( cosTh );

%% 3. Solve for parabolic mean anomaly using the cubic expression for it
Mp = .5*tan(theta/2) + (1/6)*(tan(theta/2))^3;

%% 4. Compute the elpased time from the parabolic mean anomaly
dT = Mp / (mu^2/h^3);

