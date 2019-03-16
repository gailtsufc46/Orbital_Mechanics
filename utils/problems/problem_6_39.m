%% Problem 6.39 from the textbook
% Given orbital elements, find minimum delta-v to reduce inclination to 0

%% Given 
mu = 398600.44;
a = 15000;
e = 0.5;
W = pi/4;
w = pi/6;
inc = 10*pi/180;

%{ 
Minimum delta-v to change inclination would be found by applying a plane
change burn at apogee. 

Find v at apogee, vA.

Then the delta-v for a pure rotation is:

dV = 2*vA*sin(dInc/2)

which is Eqn. 6.23 in the textbook.

%}

rA = a*(1+e);
rP = a*(1-e);
h = sqrt(2*mu)*sqrt(rP*rA/(rP+rA));
vA = h/rA;

dV = 2*vA*sin(10*pi/180/2)

