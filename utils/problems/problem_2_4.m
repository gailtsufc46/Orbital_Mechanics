% Problem 2.4
%
% At a given instant, t0, a 1000 kg Earth-orbiting satellite has the
% inertial position [3207; 5459; 2714] km and velocity [-6.532; 0.7835;
% 6.142] km/s. Solve Eqn (2.22) numerically to find the maximum altitude
% reached by the satellite and the time at which it occurs.

Re = 6378; % km                       - Earth equatorial radius
mu = 398600; % km^3/s^2               - Earth grav. parameter
r0 = [3207; 5459; 2714]; % km         - Initial position
v0 = [-6.532; 0.7835; 6.142]; % km/s  - Initial velocity

% We could make a separate m-file for this function, or we can just make an
% "anonymous" function here, without having to create/save a new file.
fun = @(t,x) [x(4:6); -mu*x(1:3)/norm(x(1:3))^3];

% Pick a good time interval. Let's use 3 hours.
T = 23*3600;
tInt = [0 T];



% Initial state vector, x0
x0 = [r0; v0];

% Compute x(t) numerically using Matlab's ODE solver
opts = odeset('abstol',1e-6,'reltol',1e-6);
[time,state] = ode45( fun, tInt, x0, opts );

% Compute the altitude
pos = state(:,1:3); % The position vector is in the first 3 columns (of every row)
vel = state(:,4:6); % The velocity vector is in the next 3 columns (of every row)
alt = sqrt(pos(:,1).^2 + pos(:,2).^2 + pos(:,3).^2)-Re;

figure, plot(time,alt)
xlabel('Time (sec)'), ylabel('Alt. (km)'), grid on, zoom on

[maxAlt,k] = max(alt);
tMaxAlt = time(k)/3600;

maxAlt, tMaxAlt



