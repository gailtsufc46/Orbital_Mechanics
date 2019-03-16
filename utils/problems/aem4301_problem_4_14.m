%{

Problem 4.14

At time t0 (relative to perigee passage), the position r0 and velocity v0  
of a satellite in the geocentric equatorial frame are
  r0 = [  -5000    -8000     -2100 ]  (km)
  v0 = [     -4        3.5      -3 ]  (km/s)

Find r and v at time t0 + 50 min.

%}

%% Method 1
%
% 1. Compute the orbital elements from position, velocity:
%     [r0,v0] -> [a,i,W,w,e,th0]
% 2. Compute the new true anomaly at DT = 50 minutes later 
%     [th0, DT] -> th1
% 3. Compute the new r,v from the orbital elements 
%     [a,i,W,w,e,th1] -> [r1,v1]

dT = 50*60; % time change in seconds

mu = 398600.4;  % EARTH'S gravitational constant (km^3/s^2)
r0 = [  -5000    -8000     -2100 ]'; % initial position (km)
v0 = [     -4;       3.5;     -3 ];  % initial velocity (km/s)

% step 1 - compute orbital elements
[a,i,W,w,e,th0] = OrbitalElementsFromRV(r0,v0,mu);

% step 2 - compute the new true anomaly 50 minutes later
E0 = 2*atan( sqrt((1-e)/(1+e)) * tan(th0/2) ); % ecc. anomaly at t0
M0 = E0 - e*sin(E0); % mean anomaly at t0
M1 = M0 + sqrt(mu/a^3)*dT; % mean anomaly at t1=t0+dT
E1 = EccAnomFromMeanAnom(M1,e,1e-6); % ecc. anomaly at t1 
% NOTE: EccAnomFromMeanAnom solves Kepler's eqn with Newton-Raphson methd.
th1 = 2*atan( sqrt((1+e)/(1-e)) * tan(E1/2) ); % true anomaly at t1

[r1,v1] = RVFromCOE(a,i,W,w,e,th1,mu); % compute [r,v] from new orbital elements


%% Method 2
% 
% 1. Compute the initial position and velocity in the perifocal frame:
%     [r0,v0] -> [r0_pq,v0_pq]

Q = GeoEqToPerifocal( i, W, w );
r0_pq = Q*r0;
v0_pq = Q*v0;

% Note: you can check this result by looking at the 3rd element of the
% above two vectors. The 3rd element ("w" direction in p-q-w frame) should
% be zero. 

% 2. Compute the initial true anomaly 
%     th0 = atan2(r0_pq(2),r0_pq(1))
th0 = atan2(r0_pq(2),r0_pq(1));

% 3. Compute the new true anomaly at DT = 50 minutes later 
%     [th0, DT] -> th1
% (Exactly the same as step 2 in Method 1)
E0 = 2*atan( sqrt((1-e)/(1+e)) * tan(th0/2) ); % ecc. anomaly at t0
M0 = E0 - e*sin(E0); % mean anomaly at t0
M1 = M0 + sqrt(mu/a^3)*dT; % mean anomaly at t1=t0+dT
E1 = EccAnomFromMeanAnom(M1,e,1e-6); % ecc. anomaly at t1 
th1 = 2*atan( sqrt((1+e)/(1-e)) * tan(E1/2) ); % true anomaly at t1

% 4. Compute the change in true anomaly
%     dTh = th1 - th0
dTh = th1 - th0;

% 5. Use Lagrange coefficients to compute new pos. and vel. in perifocal frame
%     [r0_pq,v0_pq], dth -> [r1_pq,v1_pq]
[f,g,fdot,gdot] = LagrangeCoeff( r0_pq, v0_pq, mu, dTh );

r1_pq = f*r0_pq + g*v0_pq;
v1_pq = fdot*r0_pq + gdot*v0_pq;

% 6. Compute the new pos. and vel. by rotating back to the ECI frame 
%     [r1_pq,v1_pq] -> [r1,v1]
r1b = Q'*r1_pq;
v1b = Q'*v1_pq;

%% Compare...
% We should see that r1 = r1b ... and v1 = v1b
fprintf(1,'\t r1 \t\t\t r1b\n')
disp([r1, r1b])
fprintf(1,'\n\t v1 \t\t\t v1b\n')
disp([v1, v1b])

%% Visualize...
Plot3DOrbit([a,i,W,w,e,th0],mu,1,dT)
