%% AEM 4301 - Example: Compute the variation of true anomaly across time.

%% Constants
mu = 398600.44; % grav constant for Earth [km^3/s^2]

%% Define the size and shape of the orbit
a = 15e3; % semimajor axis in [km]
e = 0.5;  % eccentricity ... use 0 < e < 1 to keep it an eccentric orbit

%% Compute the period "T" and mean orbit rate "n"
n = sqrt(mu/a^3);
T = 2*pi/n;

%% Define a time vector. t=0 corresponds to perigee crossing
t = linspace(0,T,60);

%% At each time point, compute MEAN ANOMALY "Me"
Me = n*t;

%% Compute the eccentric anomaly
f = @(E) E - e*sin(E) - Me;
df = @(E) 1 - e*cos(E);
E0 = Me;
E = zeros(size(t));
for i = 1:length(t)
  E0 = Me(i);
  f = @(x) x - e*sin(x) - Me(i);
  df = @(x) 1 - e*cos(x);
  E(i) = NewtRaph( f, df, E0, 1e-6, 100 );
end

%% Compute the true anomaly
theta = 2*atan( sqrt( (1+e)/(1-e) ) * tan(E/2) );

%% Plot angles vs time
r2d = 180/pi;
figure, plot(t,Me*r2d, t,E*r2d, t,theta*r2d, '-s', 'linewidth',2 ), grid on
legend('Mean Anomaly','Eccentric Anomaly','True Anomaly')

%% Plot path in perifocal frame
r = a*(1-e*e)./(1+e*cos(theta)); % orbit equation
x = r.*cos(theta);
y = r.*sin(theta);
figure, plot(x,y,'-o','linewidth',2), axis equal, grid on
hold on
ang = 0:.01:2*pi;
fill(Re*cos(ang),Re*sin(ang),'b')
