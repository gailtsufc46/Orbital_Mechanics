%% Problem 6.17 from the textbook
% Compute delta-v and transfer time for a bi-elliptic maneuver

Re = 6378.14;
mu = 398600.44;

e2 = 0.3;
rA = Re+300;
rC = Re+3000;
rB = rA*(1+e2)/(1-e2);

d = OrbMvrBiElliptic( rA, rB, rC, mu );

total_delta_v = d.dVTot
transfer_time = d.dTTot/3600
