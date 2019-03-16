function stateDot = RelOrbDynRHS( dr, dv, mu, h, e, th0, t )

dx = dr(1);
dy = dr(2);
dz = dr(3);

dxdot = dv(1);
dydot = dv(2);
%dzdot = dv(3);

% d/dt of relative position is just the relative velocity
rDot = dv;

% true anomaly at this time "t"
a = h*h/mu/(1-e*e);
theta = TrueAnomFromTime(t,a,e,mu,th0);

% reference orbit radius at this time
R = h^2/mu / (1+e*cos(theta) );

% d/dt of relative velocity vector:
vr = mu/h*e*sin(theta);
xDD = (2*mu/R^3 + h*h/R^4)*dx - 2*(vr*h/R^3)*dy + (2*h/R^2) * dydot;
yDD = (h*h/R^4 - mu/R^3)*dy + 2*(vr*h/R^3)*dx - (2*h/R^2) * dxdot;
zDD = -mu/R^3 * dz;

vDot = [xDD;yDD;zDD];

stateDot = [rDot; vDot];