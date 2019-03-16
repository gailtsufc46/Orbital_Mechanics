% Problem 2.5
%
% Given the initial position and velocity of an Earth-orbiting satellite,
% find its speed and distance from Earth 24 hrs later.

Re = 6378; % km             - Earth equatorial radius
mu = 398600; % km^3/s^2     - Earth grav. parameter
r0 = [6600; 0; 0]; % km     - Initial position
v0 = [0; 12; 0]; % km/s     - Initial velocity

% We could make a separate m-file for this function, or we can just make an
% "anonymous" function here, without having to create/save a new file.
fun = @(t,x) [x(4:6); -mu*x(1:3)/norm(x(1:3))^3];

% Pick a time interval of 24 hours
T = 24*3600;
tInt = [0 T];

% Initial state vector, x0
x0 = [r0; v0];

% Compute x(t) numerically using Matlab's ODE solver
[time,state] = ode45( fun, tInt, x0 );

% Compute the altitude
pos = state(:,1:3); % The position vector is in the first 3 columns (of every row)
vel = state(:,4:6); % The velocity vector is in the next 3 columns (of every row)

figure, plot(time,pos)
xlabel('Time (sec)'), ylabel('Position (km)'), grid on, zoom on

% distance from Earth center at 24 hours (last time point)
norm(pos(end,:))

% speed at 24 hours (last time point)
norm(vel(end,:))


