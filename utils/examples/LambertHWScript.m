% AEM 4301 - Lamberts Problem Homework Assignment
%

%% Case 1
clear
r1 = [8000;-3000;0];
r2 = [9000;2000;3000];
dT = 15*60;
mu = 398600.44;

[v1,v2] = LambertSolver(r1,r2,dT,mu,'pro');
el = OrbitalElementsFromRV( r1, v1, mu );
DisplayElements(el)
f = Plot3DOrbit( el, mu, 1, dT );


%% Case 2
clear
r1 = [0;1000;12000];
r2 = [500;8000;-500];
dT = 105*60;
mu = 398600.44;

[v1,v2] = LambertSolver(r1,r2,dT,mu,'pro');
el = OrbitalElementsFromRV( r1, v1, mu );
DisplayElements(el)
f = Plot3DOrbit( el, mu, 1, dT );

%% Case 3
clear
r1 = [-9313.62230551458;-3913.62230551458;3818.37661840736];
r2 = [16793.9852067806;4980.45002523689;-8353.43083665539];
dT = 18000;
mu = 398600.44;

[v1,v2] = LambertSolver(r1,r2,dT,mu,'pro');
el = OrbitalElementsFromRV( r1, v1, mu );
DisplayElements(el)
f = Plot3DOrbit( el, mu, 1, dT );

[v1,v2] = LambertSolver(r1,r2,dT,mu,'retro');
el2 = OrbitalElementsFromRV( r1, v1, mu );
DisplayElements(el2)
f = Plot3DOrbit( el2, mu, 1, dT );



%% Case 4 - Hyperbolic
clear
r1 = [23219.4878700106;8468.81579981098;3759.59332951088];
r2 = [13434.6978881845;16585.0859060679;21724.7792173688];
dT = 3600;
mu = 398600.44;

[v1,v2] = LambertSolver(r1,r2,dT,mu,'pro');
el = OrbitalElementsFromRV( r1, v1, mu );
DisplayElements(el)
f = Plot3DOrbit( el, mu, 1, dT );

%% LambertSolver
type LambertSolver

%% OrbitalElementsFromRV
type OrbitalElementsFromRV
