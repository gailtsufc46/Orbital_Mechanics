%% Problem 6.21 from the textbook
% Chase maneuver

Re = 6378.14;
mu = 398600.44;

% orbit 1
rP = 8000;
rA = 13000;
e = (rA-rP)/(rA+rP);
a = 0.5*(rA+rP);
h = sqrt(2*mu)*sqrt(rP*rA/(rP+rA));

% time from P to C
dTPC = TimeSincePeriapsis(pi/6,a,e,mu);

% time from P to D
dTPD = TimeSincePeriapsis(pi/2,a,e,mu);

% time from C to D
dTCD = dTPD-dTPC;

% radius at D
rD = h^2/mu/(1+e*cos(pi/2));

% solve Lambert problem
r1 = [rP;0;0];
r2 = [0;rD;0];
[v1T,v2T] = LambertSolver(r1,r2,dTCD,mu,'pro');

% velocity vectors at points P, D
vP = [0;h/rP;0];
vD = mu/h*[ -sin(pi/2); e+cos(pi/2); 0];

% delta-v vectors at points P, D
dVP = v1T-vP;
dVD = vD-v2T;

dVTot = norm(dVP)+norm(dVD)
