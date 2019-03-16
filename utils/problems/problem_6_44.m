%% Problem 6.44 from the textbook
% Increase altitude and change inclination using 3 different methods:
%  a) Plane change made after insertion into 600 km orbit
%  b) plane change and insertion into 600 km orbit are accomplished
%  simultaneously
%  c) plane change is made upon departing the lower orbit 

Re = 6378.14; %.14;
mu = 398600.44;

r1 = Re+300;
r2 = Re+600;
di = 20*pi/180;

v1c = sqrt(mu/r1);
v2c = sqrt(mu/r2);

%%  a) Plane change made after insertion into 600 km orbit

dh = OrbMvrHohmann( r1, r2, r1, r2, mu );
dV1 = dh.dVA;

dp = OrbMvrPlaneChangeAtApogee( dh.vBT, v2c, di, 'speed-plane' );
dV2 = dp.dVTot;

dVa = dV1+dV2;

%%  b) plane change and insertion into 600 km orbit done at same time

dh = OrbMvrHohmann( r1, r2, r1, r2, mu );
dV1 = dh.dVA;

dp = OrbMvrPlaneChangeAtApogee( dh.vBT, v2c, di, 'simultaneous' );
dV2 = dp.dVTot;

dVb = dV1+dV2;


%%  c) plane change is made upon departing the lower orbit 

dh = OrbMvrHohmann( r1, r2, r1, r2, mu );

dp = OrbMvrPlaneChangeAtApogee( v1c, dh.vAT, di, 'simultaneous' );
dV1 = dp.dVTot;
dV2 = dh.dVB;

dVc = dV1+dV2;


