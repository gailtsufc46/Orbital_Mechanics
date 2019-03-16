%% Problem 8.12 from the textbook

%% Constants
muSun = 132.7e9;
muJup = 126.7e6;
Re    = 149.6e6;
Rj    = 778.6e6;
rj    = 71490;

rp = rj+200000;

aT = .5*(Re+Rj);
Vj = sqrt(muSun/Rj);
hT = sqrt(2*muSun)*sqrt(Rj*Re/(Rj+Re));
VTj = hT/Rj;

vinf = Vj - VTj;

% flyby
e = 1+rp*vinf^2/muJup;
delta = 2*asin(1/e);
phi1 = pi;
phi2 = phi1+delta;

% new heliocentric velocity
V2rad = -vinf*sin(phi2);
V2trans = Vj+vinf*cos(phi2);
V2      = sqrt(V2rad^2+V2trans^2);

h2 = Rj*V2trans;

a2 = -muSun/2 / (V2^2/2-muSun/Rj);

e2 = sqrt(1-h2^2/a2/muSun);

%-or-%
e2sinTh2 = V2rad*h2/muSun;
e2cosTh2 = h2^2/muSun/Rj-1;
th2 = atan2( e2sinTh2, e2cosTh2 );
e2b  = e2sinTh2/sin(th2);


