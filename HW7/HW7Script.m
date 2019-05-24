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
% Define orbit paramteres
r_parking = Re+300; r_final = Re+600; v_parking = sqrt(mu_e/r_parking);
delta2 = 20*pi/180; 

% a)
% Find vt1, vt2, and v_final
ht = he(r_parking,r_final,mu_e);
vt1 = ht/r_parking;
vt2 = ht/r_final;
v_final = sqrt(mu_e/r_final);

% Find delta_v for pure rotation plane change
delta_v4 = 2*v_final*sin(delta2/2);

% Calculate total delta_v
delta_vtot = abs(vt1-v_parking)+abs(v_final-vt2)+delta_v4;
fprintf('The total delta_v for the manuever is %.4f km/s\n',delta_vtot);

% b)
% Calculate delta_v for a combined inertion and plane change
delta_v5 = sqrt((v_final-vt2)^2+4*vt2*v_final*sin(delta2/2)^2);

% Calculate total delta_v
delta_vtot2 = abs(vt1-v_parking)+delta_v5;
fprintf('The total delta_v for the manuever is %.4f km/s\n',delta_vtot2);

% c)
% Calculate delta_v for combining a plane change with the lower orbit
% departure
delta_v6 = sqrt((vt1-v_parking)^2+4*vt1*v_parking*sin(delta2/2)^2);

% Calculate total delta_v
delta_vtot3 = delta_v6+abs(v_final-vt2);
fprintf('The total delta_v for the manuever is %.4f km/s\n',delta_vtot3);
%% 5
% Orbit 1
r1circ = 8000;
v1circ = sqrt(mu_e/r1circ);

% Orbit 2
r2circ=12000; inc2 = 30*pi/180;
v2circ = sqrt(mu_e/r2circ);


% Comput the transfer orbit velocities
ht2 = he(r1circ,r2circ,mu_e);
vt21 = ht2/r1circ;
vt22 = ht2/r2circ;

% a)
% The manuever that would yield the minmum delta_v is a departure from
% orbit 1 one followed by a combination of a plane change and insertion 
% into orbit 2. This is due to the fact that the delta_v for the insertion
% is slightly smaller than the delta_v for departure, thus the delta_v 
% required to rotate the velocity vector and insert will be smaller.

% b)
% Compute delta_v for rotation and insertion into orbit 2
delta_v7 = sqrt((v2circ-vt22)^2+4*vt22*v2circ*sin(inc2/2)^2);

% Compute total delta_v for manuever
delta_vtot4 = abs(vt21-v1circ)+delta_v7;
fprintf('The total delta_v for the manuever is %.4f km/s\n',delta_vtot4);

%% 6
% a)
Lat = 44.9778; % Lattitude of Minneapolis
inc3 = 60;
Azimuth = asin(cosd(inc3)/cosd(Lat))*180/pi;
fprintf('The launch azimuth angle should be %.4f degrees\n',Azimuth);

% b) 
% An orbit with a 30 degree inclination cannot be acheived because that
% angle is less than the lattitude of Minneapolis


























