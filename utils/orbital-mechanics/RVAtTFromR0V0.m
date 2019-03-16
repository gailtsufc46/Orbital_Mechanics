function [r,v] = RVAtTFromR0V0( r0, v0, t, mu )

% Compute position and velocity at time t, given initial position and
% velocity
%
%   Inputs:
%     r0    (3,1) Initial position vector (km)
%     v0    (3,1) Initial velocity vector (km/s)
%     t     (1,n) Time elapsed since initial state (sec)
%     mu    (1,1) Gravitational constant (km^3/s^2)
%
%   Outputs:
%     r     (3,n) Position vector over time
%     v     (3,n) Velocity vector over time


% Compute "r" and "v" over a vector of time points "t"
% given the initial position and velocity, "r0" and "v0".

% compute the common orbital elements for r0, v0
[a,inc,W,w,e,th0] = OrbitalElementsFromRV( r0, v0, mu );

% compute the true anomaly values over the given time values "t"
th = TrueAnomFromTime(t,a,e,mu,th0);

dTh = th - th0;
[f,g,fd,gd] = LagrangeCoeff(r0,v0,mu,dTh);
xI = f*r0(1) + g*v0(1);
yI = f*r0(2) + g*v0(2);
zI = f*r0(3) + g*v0(3);

% compute the position and velocity from the COE at each true anomaly
%[r,v] = RVFromCOE( a,inc,W,w,e,th, mu );

r = [f*r0(1)+g*v0(1); f*r0(2)+g*v0(2); f*r0(3)+g*v0(3)];
v = [fd*r0(1)+gd*v0(1); fd*r0(2)+gd*v0(2); fd*r0(3)+gd*v0(3)];


