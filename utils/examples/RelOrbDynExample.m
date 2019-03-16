%% Numerically integrate the equations of motion for relative orbit dynamics
% Similar to Example 7.3

%% Constants
mu = 398600;
Re = 6378;

%% parameters

hP = 300; % perigee altitude (km)
e = 0.01; % eccentricity
th0 = 0;


%% reference orbit
rP = Re+hP;
a  = rP/(1-e);
rA = a*(1+e);
h = sqrt(a*mu*(1-e*e));
n = sqrt(mu/a^3);
T = 2*pi/n;

%% relative orbit initial conditions
dr0 = [-1;0;0];
dv0 = [0;2*n;0];

%% Numerically integrate the ODE
rhs = @(t,x) RelOrbDynRHS( x(1:3), x(4:6), mu, h, e, th0, t );
[tout,xout] = ode45(rhs,[0:10:(5*T)],[dr0;dv0],odeset('abstol',1e-10,'reltol',1e-10));

figure, plot(xout(:,2),xout(:,1)), grid on, axis equal
xlabel('y (Along-Track, km)')
ylabel('x (Radial, km)')
