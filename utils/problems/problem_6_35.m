%% Problem 6.35 from the textbook
% Chase maneuver with Lambert's problem
% Similar to 6.21

%% Given

mu = 398600.44;

%% orbit details
rP = 7000;
rA = 10000;
e = (rA-rP)/(rA+rP);
a = 0.5*(rA+rP);
h = sqrt(2*mu)*sqrt(rP*rA/(rP+rA));

%% Time of flight

% time from perigee to C
dTPC = TimeSincePeriapsis(120*pi/180,a,e,mu);

% time from perigee to A (half orbit period)
dTPA = OrbPeriod(a,mu)/2;

% time of flight from C to A
tOF = dTPA - dTPC;

%% Initial and Final Position Vectors

% radius at B
rB = h^2/mu/(1+e*cos(pi/2));

% solve Lambert problem
r1 = [0;rB;0];
r2 = [-rA;0;0];
[v1T,v2T] = LambertSolver(r1,r2,tOF,mu,'pro');

% velocity vectors at points C, A
vC = mu/h*[ -sin(pi/2); e+cos(pi/2); 0 ];
vA = mu/h*[ -sin(pi);   e+cos(pi);   0 ];

% delta-v vectors at points C, A
dVC = v1T-vC;
dVA = vA-v2T;

dVTot = norm(dVC)+norm(dVA)
