% HW7 Script
% Orbital Mechanics
% Garrett Ailts

%% Load Orbital HW Tool
load OrbitalHWTool

%% 6.35
% Define params
rp = 7000;
ra = 10000;
h1 = he(rp,ra,mu_e);
thetaB = pi/2;
thetaC = 120*pi/180;
thetaRend = pi;
e1 = (ra-rp)/(ra+rp);
a1 = rp/(1-e1);

% Calculate velocities at B and rendevous for orbit 1
vB = (mu_e/h1)*[-sin(thetaB); e1+cos(thetaB); 0];
vrend = (mu_e/h1)*[-sin(thetaRend); e1+cos(thetaRend); 0];

% Calculate TOF between point C and rendevous for spacecraft C
Ec = theta2E(thetaC,e1);
Mec = E2Me(Ec,e1);
n1 = sqrt(mu_e/a1^3);
TOF = (pi-Mec)/n1;

% Calculate r1 and r2 in perifocal frame for spacecraft B
r1 = (h1^2/mu_e)*(1/(1+e1*cos(thetaB)))*[cos(thetaB); sin(thetaB); 0];
r2 = (h1^2/mu_e)*(1/(1+e1*cos(pi)))*[cos(pi); sin(pi); 0];

% Solve Lambert's problem to obtain vchase1 and vchase2 for spacecraft B
[vchase1, vchase2, ~] = solveLambert(r1,r2,TOF,mu_e,'prograde');

% Obtain total delta_v for chase manuever
delta_v = norm(vrend-vchase2)+norm(vchase1-vB);
fprintf('The total delta_v for chase manuever is %.4f km/s\n',delta_v);

%% 6.39
% Define Orbit Params
a = 15e3; e = 0.5; OMEGA = 45*pi/180; omega = 30*pi/180; inc = 10*pi/180;
theta = pi; % manuever should be done at apoapse for min delta_v

% Convert orbital elements to r and v for inc=10 and inc=0 orbits
[~, v1vec, ~] = orbEl2rv(a, e, theta, OMEGA, omega, inc, mu_e);
[~, v2vec, ~] = orbEl2rv(a, e, theta, OMEGA, omega, 0, mu_e);
v1 = norm(v1vec);
v2 = norm(v2vec);

% Calculate delta_v for rotation about common apse line
delta_v1 = sqrt((v2-v1)^2+4*v1*v2*sin(inc/2)^2);
fprintf('The total delta_v for reducing the inclination to 0 is %.4f km/s\n',delta_v1);

%% 6.40
% Define Orbit Params
rc = 400+Re; inc_c = 60*pi/180;
e2 = 0.5; inc_e = 40*pi/180;

% Calculate circular orbit and elliptical orbit velocities
vc = sqrt(mu_e/rc);
ae = rc/(1-e2);
ve = sqrt(mu_e*(2/rc-1/ae));
delta = inc_c-inc_e;

delta_v3 = sqrt((ve-vc)^2+4*vc*ve*sin(delta/2)^2);
fprintf('The total delta_v for switching to the specified elliptical orbit is %.4f km/s\n',delta_v3);


%% 6.44












